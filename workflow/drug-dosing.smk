import helpers as hp

project = hp.Project("drug-dosing")


rule plot:
    input:
        project.paths.data / f"DM2470A-celltiter.txt",
        project.paths.data / f"DM2470B-celltiter.txt",
    output:
        *[
            project.paths.outs / f
            for f in (
                "combined-dose-response-curves-unconstrained.svg",
                "combined-dose-response-curves.svg",
                "dose-response-curve-231-1KB3-EN2.svg",
                "dose-response-curve-231-1KB3-EP2.svg",
                "dose-response-curve-231-1KB3.svg",
                "dose-response-curve-MDA-MB-231.svg",
                "dox-response-all-samples-curves-only-unconstrained.svg",
                "dox-response-all-samples-curves-only.svg",
                "IC50-bar.svg",
            )
        ],
    shell:
        "python drug-dosing/notebooks/response.py"
