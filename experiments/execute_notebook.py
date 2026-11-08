"""Execute every notebook cell and save its figures and tables in place."""
import os
import asyncio
from pathlib import Path

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
if os.name == 'nt':
    asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
path = ROOT / 'notebooks/STAT236_Sparse_Submatrix_Experiments.ipynb'
notebook = nbformat.read(path, as_version=4)


def on_cell_start(cell, cell_index, **kwargs):
    if cell.cell_type == 'code':
        print(f'Running cell {cell_index}: {cell.source.splitlines()[0][:100]}', flush=True)


client = NotebookClient(notebook, timeout=1800, kernel_name='python3',
                        resources={'metadata': {'path': str(ROOT)}},
                        on_cell_start=on_cell_start)
try:
    client.execute()
finally:
    nbformat.write(notebook, path)
print('Saved fully executed notebook.', flush=True)
