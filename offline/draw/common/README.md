Common layout
- io/: data access helpers (ROOT/CSV/Parquet readers, small loaders/writers, path resolution, light caching).
- plot/: plotting helpers (styles, palettes/markers, figure recipes like stack+ratio/closure/efficiency, save/export utilities).
- physics/: domain logic (selections, corrections/reweights, observable calculations, physics-driven binning, systematic helpers; keep experiment-specific here).
- utils/: generic helpers (math helpers, timing/logging, arg parsing, memoization, type/validation). Keep focused; split if it grows.
- config/: centralized configuration (constants such as luminosity/run ranges, bin edges, label maps, file lists, plot presets, YAML/JSON configs, environment switches like data/MC/tag).