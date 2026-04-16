window.BENCHMARK_DATA = {
  "lastUpdate": 1776360251131,
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
      }
    ]
  }
}