// build.rs — Compile vendored REBOUND + ASSIST C sources and our helpers.

fn main() {
    let rebound_src = "vendor/rebound/src";
    let assist_src = "vendor/assist/src";

    // REBOUND C sources (all except communication_mpi.c which needs MPI headers)
    let rebound_sources: Vec<String> = [
        "binarydiff",
        "boundary",
        "collision",
        "derivatives",
        "display",
        "fmemopen",
        "frequency_analysis",
        "glad",
        "gravity",
        "input",
        "integrator",
        "integrator_bs",
        "integrator_eos",
        "integrator_ias15",
        "integrator_janus",
        "integrator_leapfrog",
        "integrator_mercurius",
        "integrator_saba",
        "integrator_sei",
        "integrator_trace",
        "integrator_whfast",
        "integrator_whfast512",
        "output",
        "particle",
        "rebound",
        "rotations",
        "server",
        "simulationarchive",
        "tools",
        "transformations",
        "tree",
    ]
    .iter()
    .map(|name| format!("{rebound_src}/{name}.c"))
    .collect();

    // Build REBOUND as a static library
    let mut rebound_build = cc::Build::new();
    rebound_build
        .include(rebound_src)
        .define("LIBREBOUND", None)
        .define("_GNU_SOURCE", None)
        .flag_if_supported("-std=c99")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-sign-compare")
        .flag_if_supported("-Wno-unknown-pragmas")
        .flag_if_supported("-Wno-missing-field-initializers")
        .opt_level(3)
        .pic(true);

    for src in &rebound_sources {
        rebound_build.file(src);
    }
    rebound_build.compile("rebound");

    // ASSIST C sources
    let assist_sources: Vec<String> = ["assist", "forces", "spk", "ascii_ephem", "tools"]
        .iter()
        .map(|name| format!("{assist_src}/{name}.c"))
        .collect();

    let mut assist_build = cc::Build::new();
    assist_build
        .include(assist_src)
        .include(rebound_src)
        .define("LIBASSIST", None)
        .define("_GNU_SOURCE", None)
        .flag_if_supported("-std=c99")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-sign-compare")
        .flag_if_supported("-Wno-unknown-pragmas")
        .opt_level(3)
        .pic(true);

    for src in &assist_sources {
        assist_build.file(src);
    }
    assist_build.compile("assist");

    // Our C helper functions for opaque reb_simulation access
    cc::Build::new()
        .include(rebound_src)
        .include(assist_src)
        .file("src/helpers.c")
        .flag_if_supported("-std=c99")
        .opt_level(3)
        .pic(true)
        .compile("assist_rs_helpers");

    // Re-run if vendored sources change
    println!("cargo:rerun-if-changed=vendor/rebound/src");
    println!("cargo:rerun-if-changed=vendor/assist/src");
    println!("cargo:rerun-if-changed=src/helpers.c");

    // Link math library on unix
    #[cfg(unix)]
    println!("cargo:rustc-link-lib=m");
}
