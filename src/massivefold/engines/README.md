# Engines' environments

`pixi.toml` describes one environment per engine (`afmassive`, `colabfold`, `alphafold3`) and `pixi.lock` locks them:
each MassiveFold version installs exactly the locked packages (`massivefold install` runs `pixi install --locked`,
which refuses to install if `pixi.lock` does not match `pixi.toml`).

To change an engine's version:

1. edit its feature in `pixi.toml` (and `AFMASSIVE_VERSION` / `ALPHAFOLD3_VERSION` in `__init__.py` for AFmassive and
   AlphaFold3, whose `run_*.py` scripts are downloaded from the same tag),
2. update the lock with the pixi version pinned in `__init__.py` (`PIXI_VERSION`), on a machine with or without GPU:
   ```bash
   CONDA_OVERRIDE_CUDA=12 pixi lock --manifest-path src/massivefold/engines/pixi.toml
   rm -rf src/massivefold/engines/.pixi
   ```
3. commit `pixi.toml` and `pixi.lock` together.

Before a release, check that the lock matches the manifest (exits with an error otherwise):
```bash
CONDA_OVERRIDE_CUDA=12 pixi lock --check --manifest-path src/massivefold/engines/pixi.toml
```
