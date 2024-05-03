import os
import numpy as np
import pandas as pd

def update_dict(default_dict, option_dict):
    """Updates default_dict with values in option_dict"""
    updated_dict = {}
    for key in default_dict:
        updated_dict[key] = default_dict[key]
    for key in option_dict:
        updated_dict[key] = option_dict[key]
    return updated_dict


def scaling_factor(x1_old, x2_old, x1_new, x2_new):
    return (x2_new - x1_new) / (x2_old - x1_old)


def change_coordinate(x, xmin_old, xmin_new, xscaling):
    return (x - xmin_old) * xscaling + xmin_new


def change_coordinate_array(x, xmin_old, xmin_new, xscaling):
    return np.array([change_coordinate(i, xmin_old, xmin_new, xscaling) for i in x])


def add_text(xmin_old, xmin_new, xscaling, ymin_old, ymin_new, yscaling,
             text_opts):
    
    ty = (text_opts['ycoord'] - ymin_old) * yscaling + ymin_new
    tx = (text_opts['xcoord'] - xmin_old) * xscaling + xmin_new
    
    text_default = {'textcolor' : 'textColor',
                    'verticalspacing' : 10}
    
    to = {k: text_opts[k] for k in set(text_opts.keys()) - set('text')}            
    text_dict = update_dict(text_default, to)
    
    a = ''
    if type(text_opts['text']) != str:
        lc = 0
        for i in text_opts['text']:
            a += (f"\\node[text={text_dict['textcolor']},right,inner sep=0pt, "
                  f"outer sep=0pt] at ({tx}, "
                  f"{ty - lc * text_dict['verticalspacing']})"
                  f"{{{i}}};\n")
            lc += 1
    else:
        a += (f"\\node[text={text_dict['textcolor']},right,inner sep=0pt, "
              f"outer sep=0pt] at ({tx}, {ty})"
              f"{{{text_opts['text']}}};\n")
    return a


def add_line(xmin_old, xmin_new, xscaling,
             ymin_old, ymin_new, yscaling,
             line_opts):
    
    line_default = {
        'linewidth' : .6,
        'linecolor' : 'lineColor',
        'samples' : 50,
        'xmin' : 0,
        'xmax' : 1,
        'linetype' : '',
        'front' : True,
        'ymin' : None,
        'ymax' : None
    }
    
    line_dict = update_dict(line_default, line_opts)
    if line_dict['linecat'] == 'func':
        xlist = np.linspace(line_dict['xmin'], line_dict['xmax'],
                            line_dict['samples'])
        ylist = [line_dict['func'](x) for x in xlist]
    elif line_dict['linecat'] == 'data':
        xlist = line_dict['data_x']
        ylist = line_dict['data_y']
    elif line_dict['linecat'] == 'regline':
        beta = get_ols_coef(line_dict['data_y'], line_dict['data_x'])
        xlist = [line_dict['xmin'], line_dict['xmax']]
        ylist = [beta[0] + beta[1] * x for x in xlist]
    elif line_dict['linecat'] == 'vertical':
        ylist = np.linspace(line_dict['ymin'], line_dict['ymax'],
                            line_dict['samples'])
        xlist = [line_dict['func'](y) for y in ylist]
    nx = change_coordinate_array(xlist, xmin_old, xmin_new, xscaling)
    ny = change_coordinate_array(ylist, ymin_old, ymin_new, yscaling)
    if line_dict['ymin'] != None:
        ymin_l = change_coordinate(line_dict['ymin'], ymin_old, ymin_new, yscaling)
    else:
        ymin_l = -1e16
    if line_dict['ymax'] != None:
        ymax_l = change_coordinate(line_dict['ymax'], ymin_old, ymin_new, yscaling)
    else:
        ymax_l = 1e16
    a = f"\\begin{{scope}}\n"
    for j in range(len(ny)-1):
        if (ny[j] >= ymin_l) and (ny[j] <= ymax_l):
            a += (f"\t\path[draw={line_dict['linecolor']},line width="
                f"{line_dict['linewidth']}, {line_dict['linetype']}]"
                f"({nx[j]},{ny[j]}) -- ({nx[j+1]},{ny[j+1]});\n")
    a += '\end{scope}\n'
    return a, line_dict


def gen_legend(xmin_old, xmin_new, xscaling, ymin_old, ymin_new, yscaling,
               legend_opts):
    xcoord, ycoord = legend_opts['xcoord'], legend_opts['ycoord']
    ty = (ycoord - ymin_old) * yscaling + ymin_new
    tx = (xcoord - xmin_old) * xscaling + xmin_new
    
    legend_default = {'verticalspacing' : 10,
                      'textdistance' : 12,
                      'legend_marker_list' : []}
    
    leg_dict = update_dict(legend_default, legend_opts)
    
    legend_marker_default = {'shape' : 'line',
                             'linewidth' : .1,
                             'linecolor' : 'lineColor',
                             'textcolor' : 'textColor',
                             'outlinecolor' : 'dataColor',
                             'fillcolor' : 'dataColor',
                             'size' : 1,
                             'linelength' : leg_dict['textdistance'] - 2,
                             'text' : ''
    }
    
    marker_dict_list = [update_dict(legend_marker_default, i)
                        for i in legend_opts['legend_marker_list']]
    a = '\\begin{scope}\n\t'
    nr_l = 0
    for leg in marker_dict_list:
        if leg['text'] != '':
            if leg['shape'] == 'line':            
                a += (f"\path[draw={leg['linecolor']},line width="
                    f"{leg['linewidth']},line join=round,{leg['linetype']}]({tx},"
                    f"{ty - nr_l * leg_dict['verticalspacing']}"
                    f") -- ({tx + leg['linelength']},"
                    f"{ty - nr_l * leg_dict['verticalspacing']});\n")
            elif leg['shape'] == 'circle':
                a += (f"\t\path[draw={leg['outlinecolor']},"
                    f"line width={leg['linewidth']} pt,fill={leg['fillcolor']}]"
                    f"({tx + leg_dict['textdistance'] / 2}, {ty - nr_l * leg_dict['verticalspacing']}) "
                    f"{leg['shape']} ({leg['size']});\n")
            a += (f"\t\\node[text=textColor,anchor=base west,"
                f"inner sep=0pt, outer sep=0pt] at "
                f"({tx + leg_dict['textdistance']},"
                f"{ty - nr_l * leg_dict['verticalspacing'] - 2.5})"
                f"{{{leg['text']}}};\n")
            nr_l += 1
    a += '\end{scope}'
    return a


def gen_ticks(axis_loc, dmin_old, dmin_new, dscaling, axis_type, ticks_opts):
    ticks_default = {
        'labelcolor' : 'labelColor',
        'labmarg' : 5,
        'linewidth' : .6,
        'tickcolor' : 'tickColor',
        'tickmarg' : 3,
        'labscale' : 1,
        'xticks' : [],
        'yticks' : []
    }
    
    ticks_dict = update_dict(ticks_default, ticks_opts)
    if axis_type == 'x':
        ticks = ticks_dict['xticks']
    elif axis_type == 'y':
        ticks = ticks_dict['yticks']
    
    if type(ticks) == list:
        tickslist = {}
        for t in ticks:
            tickslist[t] = change_coordinate(t, dmin_old, dmin_new, dscaling)
    elif type(ticks) == dict:
        tickslist = {}
        for key, val in ticks.items():
            tickslist[val] = change_coordinate(key, dmin_old, dmin_new, dscaling)
        
    a = '\\begin{scope}\n'
    for tick, val in tickslist.items():
        vt = f"\\np{{{tick}}}"
        if axis_type == 'x':
            q1 = f"({val}, {axis_loc})"
            q2 = f"({val}, {axis_loc-ticks_dict['tickmarg']})"
            p = "below"
        elif axis_type == 'y':
            q1 = f"({axis_loc}, {val})"
            q2 = f"({axis_loc-ticks_dict['tickmarg']}, {val})"
            p = "left"
        a += (f"\t\\node[text={ticks_dict['labelcolor']},{p},inner sep="
              f"{ticks_dict['labmarg']} pt, scale={ticks_dict['labscale']}]"
              f" at {q1} {{{vt}}};\n"
              f"\t\\path[draw={ticks_dict['tickcolor']},"
              f"line width={ticks_dict['linewidth']} pt,line join=round] "
              f"{q2} -- {q1};\n")
    a += '\end{scope}\n'
    return a
    

def gen_axes(xmin, xmax, ymin, ymax, axis_opts):
    axis_default = {
        'axistype' : 'rectangle',
        'axiscolor' : 'axisColor',
        'linewidth' : .6,
        'textcolor' : 'textColor',
        'labscale' : 1,
        'xtitlerotate' : 0,
        'ytitlerotate' : 0,
        'ytitledirx' : -10,
        'xtitledirx' : 0,
        'ytitlediry' : 0,
        'xtitlediry' : -10,
        'xtitle' : '',
        'ytitle' : '',
        'xtitlepos' : 1,
        'ytitlepos' : 1
    }
    
    ax_dict = update_dict(axis_default, axis_opts)
    
    if ax_dict['axistype'] == 'rectangle':
        q = (f"\t\\path[draw={ax_dict['axiscolor']},line width="
             f"{ax_dict['linewidth']}] ({xmin}, {ymin}) rectangle ++({xmax}, {ymax});")
    else:
        if ax_dict['axistype'] == 'arrow':
            v = '-stealth'
        elif ax_dict['axistype'] == 'southwest':
            v = ''
        q = (f"\t\\path[draw={ax_dict['axiscolor']},{v},"
             f"line width={ax_dict['linewidth']}pt,line join=round]"
             f"({xmin}, {ymin}) -- ({xmax}, {ymin});\n"
             f"\t\\path[draw={ax_dict['axiscolor']},{v},"
             f"line width={ax_dict['linewidth']}pt,line join=round]"
             f"({xmin}, {ymin}) -- ({xmin}, {ymax});\n")
        
    xtitle_coord = ax_dict['xtitlepos'] * (xmax - xmin)
    ytitle_coord = ax_dict['ytitlepos'] * (ymax - ymin)
    
    a = (f"\\begin{{scope}}\n{q}"
         f"\t\\node[text={ax_dict['textcolor']},inner sep=0pt,"
         f"outer sep=0pt, scale={ax_dict['labscale']}, rotate={ax_dict['xtitlerotate']}] at "
         f"({xtitle_coord + ax_dict['xtitledirx']}, {ymin + ax_dict['xtitlediry']}){{{ax_dict['xtitle']}}};\n"
         f"\t\\node[text={ax_dict['textcolor']},inner sep=0pt,"
         f"outer sep=0pt, scale={ax_dict['labscale']}, rotate={ax_dict['ytitlerotate']}] at "
         f"({xmin + ax_dict['ytitledirx']}, {ytitle_coord + ax_dict['ytitlediry']}){{{ax_dict['ytitle']}}};\n"
         f"\end{{scope}}\n")
    
    return a


def make_line(ylims,
              xlims_tikz=(0, 300), ylims_tikz=(0, 200),
              axis_opts={}, scatter_opts={}, ticks_opts={},
              line_opts={}, text_opts={}, legend_opts={},
              shade_opts={},
              xlims=(0, 1), body_only=False):
    
    xmin_new, xmax_new = xlims_tikz
    ymin_new, ymax_new = ylims_tikz
    xmin_old, xmax_old = xlims
    ymin_old, ymax_old = ylims
        
    xscaling = scaling_factor(xmin_old, xmax_old, xmin_new, xmax_new)
    yscaling = scaling_factor(ymin_old, ymax_old, ymin_new, ymax_new)
    
    if body_only:
        a = ''
    else:
        a = ('\\begin{tikzpicture}[x=1pt,y=1pt,scale=1,'
         'framed,background rectangle/.style={draw=frameColor}]\n')
    legend_counter = 0
    if type(shade_opts) != dict:
        for s in shade_opts:
            q = shade_region(xmin_old, xmin_new, ymin_old, ymin_new,
                             xscaling, yscaling, s)
            a += q
    elif shade_opts != {}:
        q = shade_region(xmin_old, xmin_new, ymin_old, ymin_new,
                         xscaling, yscaling, shade_opts)
        a += q
    if type(line_opts) != dict:
        for i in line_opts:
            q, ldict = add_line(xmin_old, xmin_new, xscaling,
                                ymin_old, ymin_new, yscaling,
                                i)
            a += q
            if legend_opts != {}:
                legend_opts['legend_marker_list'][legend_counter] = update_dict(
                    legend_opts['legend_marker_list'][legend_counter], ldict)
                legend_counter += 1
    elif line_opts != {}:
        q, ldict = add_line(xmin_old, xmin_new, xscaling,
                            ymin_old, ymin_new, yscaling,
                            line_opts)
        a += q
        if legend_opts != {}:
            legend_opts['legend_marker_list'][legend_counter] = update_dict(
                legend_opts['legend_marker_list'][legend_counter], ldict)
            legend_counter += 1
        
    a += gen_ticks(ymin_new, xmin_old, xmin_new, xscaling, 'x', ticks_opts)
    a += gen_ticks(xmin_new, ymin_old, ymin_new, yscaling, 'y', ticks_opts)
    a += gen_axes(xmin_new, xmax_new, ymin_new, ymax_new, axis_opts)

    try:
        if scatter_opts['cmap_opts']['legend'] == 'yes':
            xmin_bar = scatter_opts['cmap_opts']['colorbardirx']
            colorbar_width = scatter_opts['cmap_opts']['colorbarwidth']
            zstart = scatter_opts['cmap_opts']['barstart']
            zmin = scatter_opts['cmap_opts']['min']
            zmax = scatter_opts['cmap_opts']['max']
            q = mpl.colormaps[scatter_opts['cmap_opts']['cmap']]
            a += gen_colorbar_legend(xmin_new + xmin_bar,
                                     xmin_new + xmin_bar + colorbar_width,
                                     zmin, zmax, zstart, q)
    except KeyError:
        pass
    if legend_opts != {}:
        a += gen_legend(xmin_old, xmin_new, xscaling,
                        ymin_old, ymin_new, yscaling, legend_opts)
    if type(text_opts) != dict:
        for i in text_opts:
            a += add_text(xmin_old, xmin_new, xscaling,
                          ymin_old, ymin_new, yscaling, i)
    elif text_opts != {}:
        a += add_text(xmin_old, xmin_new, xscaling,
                      ymin_old, ymin_new, yscaling, text_opts)
    if body_only:
        pass
    else:
        a += '\n\end{tikzpicture}'
    return a


def write_file(filename, plot):
    with open(filename, 'w', encoding='utf8') as f:
        f.write(plot)


def gen_title(xcoord, ycoord, title):
    return (f"\\begin{{scope}}\n\t\\node at ({xcoord},{ycoord})"
            f" {{{title}}};\n\\end{{scope}}\n")


def help_lines(ymin, ymax, xmin, xmax, ygrid, xgrid,
               highlight_horizontal=[]):
    horizontal_lines = [{'linecat' : 'data',
                         'data_x' : [xmin, xmax],
                         'data_y' : [l, l],
                         'linecolor' : 'labelColor',
                         'linetype' : 'densely dotted',
                         'linewidth' : .1} for l in ygrid]
    vertical_lines = [{'linecat' : 'data',
                        'data_x' : [l, l],
                        'data_y' : [ymin, ymax],
                        'linecolor' : 'labelColor',
                        'linetype' : 'densely dotted',
                        'linewidth' : .1} for l in xgrid]
    if highlight_horizontal != []:
        highlight_line = [{'linecat' : 'data',
                           'data_x' : [xmin, xmax],
                           'data_y' : [highlight_horizontal, highlight_horizontal],
                           'linecolor' : 'labelColor',
                           'linetype' : 'dashed',
                           'linewidth' : .4}]
        return horizontal_lines + vertical_lines + highlight_line
    else:
        return horizontal_lines + vertical_lines 


def gen_legend(blcoords, trcoords, lcolors, lnames, ldist=20, llength=30, lmarker=True):
    l = len(lnames)
    midpoint = (trcoords[1] - blcoords[1]) / 2
    l2 = int(np.ceil(l / 2))
    if l % 2 == 1:
        spl, spu = midpoint + blcoords[1], midpoint + blcoords[1]
    else:
        spl, spu = midpoint - ldist / 2 - + blcoords[1], midpoint + ldist / 2 - + blcoords[1]
    vals = []
    for i in range(l2):
        vals.append(spu + i * ldist)
        vals.append(spl - i * ldist)
    
    xstart_l = blcoords[0] + 10
    xend_l = blcoords[0] + 10 + llength
    vals = sorted(set(vals), reverse=True)
    t = f"\\begin{{scope}}\n"
    for i in range(l):
        t += (f"\t\\path[draw={lcolors[i]},line width=1, solid]"
              f"({xstart_l},{vals[i]}) -- ({xend_l},{vals[i]});\n")
        t += f"\t\\node[text=labelColor,right,inner sep=5 pt, scale=1] at ({xend_l},{vals[i]}) {{{lnames[i]}}};\n"
    if lmarker == True:
        t += f"\t\\node[fill=lcolor1,circle,inner sep=1.5pt] at ({(xstart_l + xend_l) / 2},{vals[0]}) {{}}\n;"
        t += f"\t\\node[fill=lcolor2,regular polygon, regular polygon sides=4,inner sep=1.5pt] at ({(xstart_l + xend_l) / 2},{vals[1]}) {{}}\n;"
        t += f"\t\\node[fill=lcolor3,regular polygon, regular polygon sides=3,inner sep=1.5pt] at ({(xstart_l + xend_l) / 2},{vals[2]}) {{}}\n;"
    t += f"\\end{{scope}}\n"
    return t

def gen_legend_horizontal(mx, my, lcolors, lnames, hdist=[], llength=30, lmarker=True):
    l = len(lnames)
    midpoint = my - 30
    if hdist == []:
        hdist = [100 for i in range(l)]
    t = f"\\begin{{scope}}\n"
    for i in range(l):
        stval = mx + 10 + i * hdist[i]
        endval = stval + llength
        t += (f"\t\\path[draw={lcolors[i]},line width=1, solid]"
              f"({stval},{midpoint}) -- ({endval},{midpoint});\n")
        t += f"\t\\node[text=labelColor,right,inner sep=5 pt, scale=1] at ({endval},{midpoint}) {{{lnames[i]}}};\n"
    if lmarker == True:
        shapelist = ['circle',
                     'regular polygon, regular polygon sides=4',
                     'regular polygon, regular polygon sides=3',
                     'star']
        for i in range(l):
            stval = mx + 10 + i * hdist[i]
            endval = stval + llength
            t += f"\t\\node[fill=lcolor{i+1},{shapelist[i]},inner sep=1.5pt] at ({(stval + endval) / 2},{midpoint}) {{}}\n;"
    t += f"\\end{{scope}}\n"
    return t


def gen_plot(xmin, xmax, ymin, ymax, nr_rows, nr_cols, Yl, Xl, lcolorlist,
             xgrid, ygrid, xticks, yticks, legend='vertical',
             highlight_horizontal=[], lmarker=[], titlelist=[],
             graph_coords = [130, 130], tl_rel_coord = [120, 100],
             ltitle=False):
    xmin_old = xmin
    xmax_old = xmax
    ymin_old = ymin
    ymax_old = ymax
    rows, cols = nr_rows, nr_cols
    blcoordlist = []
    nry = Yl[0].shape[1]

    for i in range(rows):
        for j in range(cols):
            x = graph_coords[0] * j
            y = -graph_coords[1] * i
            blcoordlist.append([x, y])

    Q = ''
    for i in range(len(Yl)):
        Y = Yl[i]
        X = Xl[i]
        DATA_LINES = []
        ym = []
        for l in range(nry):
            y = Y[:, l]
            data_lines = {'data_x' : X,
                          'data_y' : y,
                          'linecat' : 'data',
                          'xmin' : xmin_old,
                          'xmax' : xmax_old,
                          'linecolor' : lcolorlist[l],
                          'linewidth' : 1}
            if lmarker != []:
                if lmarker[i] != []:
                    if lmarker[i][l] == []:
                        yval = np.nan
                    else:
                        yval = y[X == lmarker[i][l]][0]
                    ym.append(yval)
            DATA_LINES.append(data_lines)

        xmin_new = blcoordlist[i][0]
        ymin_new = blcoordlist[i][1]
        xmax_new = xmin_new + tl_rel_coord[0]
        ymax_new = ymin_new + tl_rel_coord[1]

        p = ''
        if lmarker != []:
            if lmarker[i] != []:
                p = f"\\begin{{scope}}\n"
                shapelist = ['circle',
                     'regular polygon, regular polygon sides=4',
                     'regular polygon, regular polygon sides=3',
                     'star']
                for l in range(len(lmarker[i])):
                    if lmarker[i][l] != []:
                        xscaling = scaling_factor(
                            xmin_old, xmax_old, xmin_new, xmax_new)
                        yscaling = scaling_factor(
                            ymin_old, ymax_old, ymin_new, ymax_new)

                        xc = change_coordinate_array(
                                [lmarker[i][l]], xmin_old, xmin_new, xscaling)
                        yc = change_coordinate_array(
                                [ym[l]], ymin_old, ymin_new, yscaling)
                        p += f"\t\\node[fill=lcolor{l+1},{shapelist[l]},inner sep=1.5pt] at ({xc[0]},{yc[0]}) {{}}\n;"
                p += f"\\end{{scope}}\n"

        hlines = help_lines(ymin_old, ymax_old, xmin_old, xmax_old, ygrid, xgrid,
                                highlight_horizontal=highlight_horizontal)
        DATA_LINES = DATA_LINES + hlines
        if i % nr_cols == 0:
            yt = yticks
        else:
            yt = []
        if i >= len(Yl) - nr_cols:
            xt = xticks
        else:
            xt = []
        q = make_line(
            ylims=(ymin_old, ymax_old), xlims=(xmin_old, xmax_old),
            xlims_tikz=(xmin_new, xmax_new), ylims_tikz=(ymin_new, ymax_new),
            axis_opts={'axistype' : 'southwest'},
            line_opts=DATA_LINES,
            body_only=True,
            ticks_opts={'xticks' : xt, 'yticks' : yt})
        if titlelist != []:
            if titlelist[i] != '':
                t = gen_title(blcoordlist[i][0] + tl_rel_coord[0] / 2,
                    blcoordlist[i][1] + tl_rel_coord[1] + 10, titlelist[i])
            else:
                t = ''
        else:
            t = ''
        Q += q
        Q += p
        Q += t
    Q = '\\begin{tikzpicture}[x=1pt,y=1pt,scale=1,framed,background rectangle/.style={draw=frameColor}]\n' + Q
    if legend == 'vertical':
        bl = blcoordlist[len(Yl)]
        tr = [bl[0] + tl_rel_coord[0], bl[1] + tl_rel_coord[1]]
        Q += gen_legend(bl, tr, lcolorlist, ['Naive', 'Conservative', 'Always valid'])
    elif legend=='horizontal':
        mx = min([i[0] for i in blcoordlist])
        my = min([i[1] for i in blcoordlist])
        Q += gen_legend_horizontal(mx, my, lcolorlist,
            ['Naive', 'Conservative', 'Always valid', 'GST'],
             hdist=[0, 75, 95, 100], lmarker=True)
    if ltitle:
        y = 120
        Ql = gen_title(-40, 80, 'DGP')
        for i in range(nr_rows):
            y -= 80
            q = gen_title(-40, y, f'{i+3}')
            Ql += q
        Q += Ql    
    Q += '\\end{tikzpicture}'
    return Q


def gen_YX_list(DF, ynames, xname, bias=False):
    Yl = []
    Xl = []
    ATE = {1 : 0, 2 : 0, 3 : .2, 4 : .2, 5 : .2, 6 : .2, 7 : .2, 8 : .2}
    for i in DF['distrYid'].unique():
        df = DF[DF['distrYid'] == i]
        if bias:
            Y = (df[ynames] - ATE[i]).to_numpy()
        else:
            Y = df[ynames].to_numpy()
        X = df[xname].to_numpy()
        Yl.append(Y)
        Xl.append(X)
    return Yl, Xl


def gen_figures(str_val):
    mFILE = f"simulations\\out_FWCID_{str_val}_agg.csv"
    DF = pd.read_csv(mFILE)
    lcolorlist = ['lcolor1', 'lcolor2', 'lcolor3']

    xmin = 0
    xmax = .52
    xgrid = [.25, .5]
    xticks = [0, .25, .5]
    nr_rows = 3
    nr_cols = 3

    lmlist = [[.25, .25, .25]]*8
    lmlist2 = [[], [], [], [], [.25, .25, .25], [], [.25, .25, .25], [.4, .4, .4]]

    for j in ['naive', 'cons', 'av']:
        DF[f'n_ratio_{j}'] = DF[f'n_{j}_mean'] / DF['nbar']

    tlist = [f"DGP {i}" for i in range(1, len(DF['distrYid'].unique()) + 1)]

    Yl, Xl = gen_YX_list(DF,
                        ['n_ratio_naive', 'n_ratio_cons', 'n_ratio_av'],
                        'd')
    Q1 = gen_plot(xmin, xmax, 0, 8.2, nr_rows, nr_cols, Yl, Xl, lcolorlist,
                xgrid, [2, 4, 6], xticks, [0, 2, 4, 6, 8], legend='vertical',
                highlight_horizontal=1, lmarker=lmlist, titlelist=tlist)

    Yl, Xl = gen_YX_list(DF,
                        ['te_naive_mean', 'te_cons_mean', 'te_av_mean'],
                        'd', bias=True)
    Q2 = gen_plot(xmin, xmax, -.07, .07, nr_rows, nr_cols, Yl, Xl, lcolorlist,
                xgrid, [-.05, 0, .05], xticks, [-.05, 0, .05], legend='vertical',
                highlight_horizontal=0, lmarker=lmlist2, titlelist=tlist)

    Yl, Xl = gen_YX_list(DF,
                        ['cov_naive_mean', 'cov_cons_mean', 'cov_av_mean'],
                        'd')
    Q3 = gen_plot(xmin, xmax, .78, 1, nr_rows, nr_cols, Yl, Xl, lcolorlist,
                xgrid, [.8, .9, 1], xticks, [.8, .9, 1], legend='vertical',
                highlight_horizontal=.9, lmarker=lmlist, titlelist=tlist)


    write_file(f'figures\\nratio_{str_val}.tex', Q1)
    write_file(f'figures\\bias_{str_val}.tex', Q2)
    write_file(f'figures\\coverage_{str_val}.tex', Q3)

    tFILE = f"simulations\\out_FPD_{str_val}.csv"
    DF = pd.read_csv(tFILE)
    DF = DF[DF['distrYid'] >=3]
    DF['ratio'] = (DF.groupby(['distrYid', 'tau_H1'])['maxn'].rank() - 1) / 2 + 1

    for i in range(4):
        DF[f'n_ratio{i+1}'] = DF[f'n{i+1}'] / DF['maxn']

    YL = []
    XL = []
    YL2 = []
    XL2 = []
    YL3 = []
    XL3 = []
    for nrat in [1, 1.5, 2]:
        DFa = DF[(DF['ratio'] == nrat) ]
        
        Yla, Xla = gen_YX_list(DFa,
            ['est1', 'est2', 'est3', 'est4'], 'tau_H1')
        Ylb, Xlb = gen_YX_list(DFa,
            ['n_ratio1', 'n_ratio2', 'n_ratio3', 'n_ratio4'], 'tau_H1')
        Ylc, Xlc = gen_YX_list(DFa,
            ['power1', 'power2', 'power3', 'power4'], 'tau_H1')
        YL.append(Yla)
        XL.append(Xla)
        YL2.append(Ylb)
        XL2.append(Xlb)
        YL3.append(Ylc)
        XL3.append(Xlc)
    YLIST = []
    XLIST = []
    YLIST2 = []
    XLIST2 = []
    YLIST3 = []
    XLIST3 = []
    for i in range(len(DF['distrYid'].unique())):
        for j in range(3):
            YLIST.append(YL[j][i] - .2)
            XLIST.append(XL[j][i])
            YLIST2.append(YL2[j][i])
            XLIST2.append(XL2[j][i])
            YLIST3.append(YL3[j][i])
            XLIST3.append(XL3[j][i])

    lcolorlist = ['lcolor1', 'lcolor2', 'lcolor3', 'lcolor4']


    xmin = .07
    xmax = .43
    xgrid = [.1, .2, .3, .4]
    xticks = [.1, .2, .3, .4]
    nr_rows = 6
    nr_cols = 3

    mlist1 = [[.22, .28, .25, .25], [.22, .28, .25, .25], [.22, .28, .2, .25],
            [.22, .28, .15, .25], [.22, .28, .15, .25], [.22, .28, .25, .25],
            [.3, .35, .15, .25], [.35, .35, .15, .25], [.35, .35, .15, .25],
            [.27, .32, .15, .25], [.27, .32, .15, .25], [.27, .32, .15, .25],
            [.39, .35, .15, .25], [.39, .35, .15, .25], [.39, .35, .15, .25],
            [.27, .32, .15, .25], [.27, .32, .15, .25], [.27, .32, .15, .25]]

    mlist2 = [[.3, .3, .25, .25], [.15, .15, .25, .28], [.3, .3, .25, .27],
            [.2, .32, .15, .25], [.22, .28, .15, .25], [.22, .22, .25, .25],
            [.3, .3, .12, .18], [.25, .25, .25, .25], [.25, .35, .25, .21],
            [.25, .35, .12, .25], [.17, .25, .2, .3], [.27, .2, .2, .15],
            [.25, .38, .12, .25], [.3, .25, .2, .2], [.25, .35, .2, .15],
            [.25, .35, .12, .25], [.3, .25, .2, .17], [.25, .35, .2, .15]]

    mlist3 = [[.15, .25, .15, .25], [.25, .3, .25, .25], [.35, .3, .2, .25],
            [.25, .25, .17, .25], [.25, .25, .25, .25], [.35, .3, .25, .25],
            [.25, .15, .17, .25], [.25, .25, .25, .25], [.35, .3, .25, .25],
            [.25, .12, .15, .25], [.25, .25, .25, .25], [.37, .39, .25, .25],
            [.25, .15, .15, .25], [.35, .32, .25, .25], [.22, .18, .25, .25],
            [.25, .15, .15, .25], [.32, .32, .25, .25], [.22, .18, .25, .25]]
    lmlist = [[.25, .25, .25, .25]]*18

    tlist = ['$N^{max}=\\bar{N}$', '$N^{max}=1.5\\bar{N}$', '$N^{max}=2\\bar{N}$',
            '', '', '', '', '', '', '', '', '', '', '', '', '', '', '']

    Q4 = gen_plot(xmin, xmax, -.06, .17, nr_rows, nr_cols, YLIST, XLIST, lcolorlist,
                xgrid, [0, .1], xticks, [0, .1], legend='horizontal',
                highlight_horizontal=0, titlelist=tlist, lmarker=mlist1,
                graph_coords = [130, 80], tl_rel_coord = [120, 70], ltitle=True)

    Q5 = gen_plot(xmin, xmax, 0, 1.02, nr_rows, nr_cols, YLIST2, XLIST2, lcolorlist,
                xgrid, [0, .5, 1], xticks, [0, .5, 1], legend='horizontal',
                titlelist=tlist, lmarker=mlist2,
                graph_coords = [130, 80], tl_rel_coord = [120, 70], ltitle=True)

    Q6 = gen_plot(xmin, xmax, 0, 1.02, nr_rows, nr_cols, YLIST3, XLIST3, lcolorlist,
                xgrid, [0, .5, 1], xticks, [0, .5, 1], legend='horizontal',
                titlelist=tlist, lmarker=mlist3,
                graph_coords = [130, 80], tl_rel_coord = [120, 70], ltitle=True)

    write_file(f'figures\\bias_test_{str_val}.tex', Q4 )
    write_file(f'figures\\nratio_test_{str_val}.tex', Q5)
    write_file(f'figures\\power_test_{str_val}.tex', Q6)


def main():
    gen_figures(args.speed)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate figures')
    parser.add_argument('speed')
    args = parser.parse_args()
    main()
