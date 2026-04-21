window.BENCHMARK_DATA = {
  "lastUpdate": 1776796709173,
  "repoUrl": "https://github.com/B612-Asteroid-Institute/assist-rs",
  "entries": {
    "assist-rs Benchmarks": [
      {
        "commit": {
          "author": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "committer": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "distinct": true,
          "id": "4fcfdc0010143751cd1dd44135275f0b57c9a803",
          "message": "Apply rustfmt to recent commits\n\nRustfmt drift across the 8 preceding commits; fix in one sweep rather\nthan amending pushed history. No functional changes.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-04-16T10:18:21-07:00",
          "tree_id": "9b08caeed1e62be43ed948b54a04c0dac182c288",
          "url": "https://github.com/B612-Asteroid-Institute/assist-rs/commit/4fcfdc0010143751cd1dd44135275f0b57c9a803"
        },
        "date": 1776360250346,
        "tool": "cargo",
        "benches": [
          {
            "name": "propagate_single/rust_api/1",
            "value": 497090,
            "range": "± 13784",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/10",
            "value": 1321596,
            "range": "± 40310",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/100",
            "value": 10155775,
            "range": "± 31514",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/without_stm",
            "value": 500329,
            "range": "± 2470",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/with_stm",
            "value": 1107083,
            "range": "± 22386",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/gravity_only",
            "value": 499343,
            "range": "± 20592",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/with_a2",
            "value": 546701,
            "range": "± 2471",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/rust_api",
            "value": 503355,
            "range": "± 6437",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/raw_c_ffi",
            "value": 501211,
            "range": "± 15667",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/serial_28_orbits",
            "value": 14200840,
            "range": "± 69861",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/rayon_28_orbits",
            "value": 10647135,
            "range": "± 67357",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/1",
            "value": 350027,
            "range": "± 1446",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/10",
            "value": 442749,
            "range": "± 3026",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/30",
            "value": 500827,
            "range": "± 1536",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/100",
            "value": 745503,
            "range": "± 13359",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/365",
            "value": 1453048,
            "range": "± 29414",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "committer": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "distinct": true,
          "id": "5d2eaf8d39f08852ca813f9b4c6d016d2c2a48f7",
          "message": "Add covariance propagation methods on PropagatedState\n\nTwo methods for linearly propagating an initial covariance to the\npropagation epoch using the STM (and, optionally, the non-grav\npartials):\n\n- propagate_covariance(&P₀_6x6)              → P(t)_6x6   = Φ·P₀·Φᵀ\n- propagate_covariance_with_nongrav(&P₀_9x9) → P(t)_6x6   = J·P₀·Jᵀ\n  where J = [Φ | G] is 6×9 and P₀ is over (x, y, z, vx, vy, vz, A1,\n  A2, A3).\n\nBoth return None when the required partials are absent (e.g., method\ncalled on a state propagated without `compute_stm`, or the non-grav\nvariant called on a gravity-only orbit). No change to the existing\nSTM output — covariance is a convenience on top, not a replacement.\n\nTests cover: identity-P₀ → Φ·Φᵀ, zero-STM path, state-only 9×9 reduces\nto 6×6 path, pure-nongrav 9×9 equals G·P_AA·Gᵀ, and full-identity 9×9\nis the sum of state and nongrav contributions.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-04-16T10:19:34-07:00",
          "tree_id": "9b08caeed1e62be43ed948b54a04c0dac182c288",
          "url": "https://github.com/B612-Asteroid-Institute/assist-rs/commit/5d2eaf8d39f08852ca813f9b4c6d016d2c2a48f7"
        },
        "date": 1776360330172,
        "tool": "cargo",
        "benches": [
          {
            "name": "propagate_single/rust_api/1",
            "value": 497585,
            "range": "± 3432",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/10",
            "value": 1325070,
            "range": "± 28830",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/100",
            "value": 10168444,
            "range": "± 67401",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/without_stm",
            "value": 501946,
            "range": "± 4567",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/with_stm",
            "value": 1108308,
            "range": "± 34889",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/gravity_only",
            "value": 502502,
            "range": "± 2683",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/with_a2",
            "value": 549230,
            "range": "± 2504",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/rust_api",
            "value": 501761,
            "range": "± 4468",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/raw_c_ffi",
            "value": 500053,
            "range": "± 1874",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/serial_28_orbits",
            "value": 14178168,
            "range": "± 55110",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/rayon_28_orbits",
            "value": 10676382,
            "range": "± 168310",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/1",
            "value": 352696,
            "range": "± 34741",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/10",
            "value": 450540,
            "range": "± 10364",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/30",
            "value": 501057,
            "range": "± 6445",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/100",
            "value": 743223,
            "range": "± 9442",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/365",
            "value": 1444493,
            "range": "± 42820",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "committer": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "distinct": true,
          "id": "7dd7336388f3b64a56d9f8ff39602c083bd56132",
          "message": "CI: use B612 data packages for observatory codes and EOP kernels\n\nThe Test job already installs naif-de440 and jpl-small-bodies-de441-n16\nfrom B612's PyPI, but the Horizons v2 test falls back to downloading\nMPC obscodes and NAIF EOP kernels over HTTP at test time. That's slow,\nnetwork-dependent, and couples CI correctness to upstream uptime.\n\nInstall the remaining B612-published data packages and point the test\nloaders at them via env vars (MPC_OBSCODES_PATH, ASSIST_EOP_*):\n\n  - mpc-obscodes           → obscodes_extended.json\n  - naif-eop-high-prec     → earth_latest_high_prec.bpc\n  - naif-eop-historical    → earth_620120_240827.bpc\n  - naif-eop-predict       → earth_200101_990827_predict.bpc\n\nThe test-side loader (added in the previous commit) already honors\nthese env vars and falls back to ASSIST_DATA_DIR for local dev, so\nnothing else needs to change.\n\nBenchmark job doesn't need the observatory or EOP kernels — propagation\nbenches only touch planetary + asteroid ephemerides — so it's left\nalone.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-04-16T10:57:56-07:00",
          "tree_id": "baaa81f726daddfcaa3a6ed7232bd9e3bec9e87f",
          "url": "https://github.com/B612-Asteroid-Institute/assist-rs/commit/7dd7336388f3b64a56d9f8ff39602c083bd56132"
        },
        "date": 1776362630786,
        "tool": "cargo",
        "benches": [
          {
            "name": "propagate_single/rust_api/1",
            "value": 521846,
            "range": "± 16423",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/10",
            "value": 1376536,
            "range": "± 27055",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/100",
            "value": 10556608,
            "range": "± 288570",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/without_stm",
            "value": 522696,
            "range": "± 2208",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/with_stm",
            "value": 1161451,
            "range": "± 28737",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/gravity_only",
            "value": 524598,
            "range": "± 1703",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/with_a2",
            "value": 574005,
            "range": "± 2749",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/rust_api",
            "value": 527440,
            "range": "± 1980",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/raw_c_ffi",
            "value": 526618,
            "range": "± 2155",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/serial_28_orbits",
            "value": 14898516,
            "range": "± 360922",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/rayon_28_orbits",
            "value": 11436529,
            "range": "± 64921",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/1",
            "value": 370756,
            "range": "± 2346",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/10",
            "value": 468775,
            "range": "± 3557",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/30",
            "value": 526475,
            "range": "± 12630",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/100",
            "value": 783665,
            "range": "± 22737",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/365",
            "value": 1521234,
            "range": "± 35369",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "committer": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "distinct": true,
          "id": "0268b50756c6edca27191c6121c0a4d216736722",
          "message": "Add Earth orientation kernels to DataManager catalog\n\nFlagged during code review but skipped when the conversation diverged\ninto the SPICE precedence work — DataManager only downloaded\nde440.bsp, sb441-n16.bsp, and obscodes_extended.json. A user relying\non the default-feature data manager had no way to obtain the EOP\nkernels and fell back to the ~50 mas GMST approximation for ground\nobservatories.\n\nRegister the three NAIF Earth orientation PCKs in DEFAULT_KERNELS:\n\n  - earth_latest_high_prec.bpc         (dynamic; weekly updates)\n  - earth_620120_250826.bpc            (historical, 1962 → 2025)\n  - earth_2025_250826_2125_predict.bpc (long-term predict → 2125)\n\n`AssistDataPaths` gains three new fields plus an `eop_kernels()`\nhelper that returns them in SPICE-idiomatic load order\n(predict → historical → current) so the high-precision kernel wins\nat epochs it covers. README example updated to show the full flow.\n\nTwo new tests enforce invariants the human eye missed the first time:\nevery kernel in `DEFAULT_KERNELS` must be reachable through\n`AssistDataPaths` (and vice versa), and `eop_kernels()` must return\nthem in last-in-wins order.\n\nNAIF periodically republishes the historical and predict kernels with\nnew date-range suffixes; when that happens, bump the two filenames in\n`DEFAULT_KERNELS` and in `paths()`. `earth_latest_high_prec.bpc` is\nNAIF's stable endpoint for continuous updates.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-04-16T23:05:38-07:00",
          "tree_id": "5e1bdc63e31047a6c01732b9ea540e8c596cf337",
          "url": "https://github.com/B612-Asteroid-Institute/assist-rs/commit/0268b50756c6edca27191c6121c0a4d216736722"
        },
        "date": 1776406348140,
        "tool": "cargo",
        "benches": [
          {
            "name": "propagate_single/rust_api/1",
            "value": 498455,
            "range": "± 50607",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/10",
            "value": 1321603,
            "range": "± 106468",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/100",
            "value": 10157060,
            "range": "± 44045",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/without_stm",
            "value": 500615,
            "range": "± 2928",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/with_stm",
            "value": 1114049,
            "range": "± 29696",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/gravity_only",
            "value": 499189,
            "range": "± 2566",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/with_a2",
            "value": 545602,
            "range": "± 3248",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/rust_api",
            "value": 501270,
            "range": "± 7440",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/raw_c_ffi",
            "value": 501510,
            "range": "± 14373",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/serial_28_orbits",
            "value": 14216164,
            "range": "± 614548",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/rayon_28_orbits",
            "value": 10703107,
            "range": "± 187914",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/1",
            "value": 349567,
            "range": "± 23716",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/10",
            "value": 442717,
            "range": "± 19848",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/30",
            "value": 502581,
            "range": "± 19647",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/100",
            "value": 742741,
            "range": "± 12641",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/365",
            "value": 1445526,
            "range": "± 32040",
            "unit": "ns/iter"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "moeyensj@gmail.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "committer": {
            "email": "moeyensj@users.noreply.github.com",
            "name": "Joachim Moeyens",
            "username": "moeyensj"
          },
          "distinct": true,
          "id": "260bcbd3610ab972e41b6ea165e98187127ee9ff",
          "message": "Bundle ephemeris + observatory table into AssistData\n\nThe five public high-level entry points (`assist_propagate`,\n`assist_propagate_single`, `assist_get_state`,\n`assist_generate_ephemeris`, `assist_generate_ephemeris_single`) and\n`PropagatorPool::new` now take a single `&AssistData` argument instead\nof an `ephem: &Ephemeris` plus (for the ephemeris-generating paths) a\nseparate `obs_table: Option<&ObservatoryTable>`.\n\n`AssistData` is a small bundle:\n\n```rust\nlet data = AssistData::new(ephem);                    // propagation only\nlet data = AssistData::new(ephem).with_observatory(obs_table);\n```\n\nCallers load these resources once at startup anyway; plumbing them\nthrough as separate args at every call was noise. The observatory\ntable can also carry an EarthOrientation via its own builder, so the\nwhole \"data dependencies\" picture now lives in one type.\n\nNo behavioral change — all tests pass bitwise-identically. README\nexamples, benches, and validation tests updated to match.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>",
          "timestamp": "2026-04-21T11:28:06-07:00",
          "tree_id": "0c2f83561f63570b6ff2b471de647d4c58d8a6f8",
          "url": "https://github.com/B612-Asteroid-Institute/assist-rs/commit/260bcbd3610ab972e41b6ea165e98187127ee9ff"
        },
        "date": 1776796708103,
        "tool": "cargo",
        "benches": [
          {
            "name": "propagate_single/rust_api/1",
            "value": 383880,
            "range": "± 1504",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/10",
            "value": 1013172,
            "range": "± 59054",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_single/rust_api/100",
            "value": 7783590,
            "range": "± 157889",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/without_stm",
            "value": 386963,
            "range": "± 2784",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_stm/with_stm",
            "value": 908688,
            "range": "± 3900",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/gravity_only",
            "value": 389205,
            "range": "± 3081",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_nongrav/with_a2",
            "value": 418264,
            "range": "± 1449",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/rust_api",
            "value": 396148,
            "range": "± 19663",
            "unit": "ns/iter"
          },
          {
            "name": "rust_vs_raw_c/raw_c_ffi",
            "value": 387306,
            "range": "± 1521",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/serial_28_orbits",
            "value": 11051604,
            "range": "± 238992",
            "unit": "ns/iter"
          },
          {
            "name": "parallel/rayon_28_orbits",
            "value": 9399545,
            "range": "± 38109",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/1",
            "value": 278018,
            "range": "± 1033",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/10",
            "value": 346585,
            "range": "± 5167",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/30",
            "value": 389725,
            "range": "± 2351",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/100",
            "value": 577485,
            "range": "± 24170",
            "unit": "ns/iter"
          },
          {
            "name": "duration_scaling/days/365",
            "value": 1115118,
            "range": "± 27232",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_30d/unpooled",
            "value": 51485458,
            "range": "± 259173",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_30d/pooled",
            "value": 51114829,
            "range": "± 390445",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_30d/unpooled_with_stm",
            "value": 117737956,
            "range": "± 1080161",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_30d/pooled_with_stm",
            "value": 117431866,
            "range": "± 371229",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_365d/unpooled",
            "value": 142167424,
            "range": "± 335118",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_365d/pooled",
            "value": 141652966,
            "range": "± 1002762",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_365d/unpooled_with_stm",
            "value": 320269507,
            "range": "± 1023897",
            "unit": "ns/iter"
          },
          {
            "name": "pool_vs_unpooled_365d/pooled_with_stm",
            "value": 320173240,
            "range": "± 1844074",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_batch/serial_loop_128",
            "value": 51009567,
            "range": "± 360660",
            "unit": "ns/iter"
          },
          {
            "name": "propagate_batch/batch_api_128",
            "value": 42962826,
            "range": "± 179954",
            "unit": "ns/iter"
          },
          {
            "name": "generate_ephemeris/earth_7_observers_30d",
            "value": 2426031,
            "range": "± 12596",
            "unit": "ns/iter"
          }
        ]
      }
    ]
  }
}