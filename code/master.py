import subprocess


def run_all(replicate_paper):
    if replicate_paper:
        subprocess.run(["julia", "code\\FWCID_fast_sim.jl", "1000", ".01", ".5", ".01"])
        #subprocess.run(["julia", "code\\gst_z.jl", ".1", ".4", ".01"])
        subprocess.run(["julia", "code\\FPD_fast_sim.jl", "1000", ".1", ".4", ".01"])
        subprocess.run(["python", "code\\gen_tikz_figures.py", "fast"])
        for figname_stem in ["bias", "coverage", "nratio",
                            "bias_test", "power_test", "nratio_test"]:
            subprocess.run(["python", "code\\convert_tikz.py",
                           figname_stem + '_fast'])
    else:
        subprocess.run(["julia", "code\\FWCID_fast_sim.jl", "10", ".3", ".5", ".01"])
        subprocess.run(["julia", "code\\FWCID_slow_sim.jl", "10", ".3", ".5", ".01"])
        subprocess.run(["julia", "code\\gst_z.jl", ".37", ".4", ".01"])
        subprocess.run(["julia", "code\\FPD_fast_sim.jl", "10", ".37", ".4", ".01"])
        subprocess.run(["julia", "code\\FPD_slow_sim.jl", "10", ".37", ".4", ".01"])


def main():
    run_all(args.replicate_paper)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Run master')
    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('-r', '--repl', dest='replicate_paper', action='store_true')
    feature_parser.add_argument('-c', '--comp', dest='replicate_paper', action='store_false')
    parser.set_defaults(replicate_paper=True)
    args = parser.parse_args()
    main()
