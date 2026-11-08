"""Export a reading copy: retain outputs and hide cells tagged hide-input."""
from pathlib import Path

import nbformat
from nbconvert import HTMLExporter
from traitlets.config import Config

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / 'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb'
config = Config()
config.TagRemovePreprocessor.remove_input_tags = {'hide-input'}
config.HTMLExporter.exclude_input_prompt = True
config.HTMLExporter.exclude_output_prompt = True
exporter = HTMLExporter(config=config)
body, _ = exporter.from_notebook_node(nbformat.read(path, as_version=4))
target = path.with_suffix('.html')
target.write_text(body, encoding='utf-8')
print('Saved reading copy:', target)
