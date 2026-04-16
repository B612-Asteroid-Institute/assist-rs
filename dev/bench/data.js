window.BENCHMARK_DATA = {
  "lastUpdate": 1776362631507,
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
      }
    ]
  }
}