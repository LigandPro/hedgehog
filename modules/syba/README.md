# SYBA (vendored)

**SYnthetic BAyesian classifier** — fragment-based score for easy- vs hard-to-synthesize
molecules (ECFP4-like fragments + learned fragment weights).

In Hedgehog it shows up as the `syba` entry in synthesis `enabled_scores`:

```yaml
enabled_scores:
  - sa
  - syba
  - rascore
syba_score_min: 0
syba_score_max: inf
```

Positive SYBA → easier; negative → harder (see the papers below).

This checkout is for local/offline use. Training notebooks are **not** vendored;
use the upstream project for those.

- Code / notebooks: [lich-uct/syba](https://github.com/lich-uct/syba)  
- Paper: [SYBA](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00439-2)  
- Related: [Nonpher](http://dx.doi.org/10.1186/s13321-017-0206-2)
