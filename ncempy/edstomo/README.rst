ncempy.edstomo
==========

The ``ncempy.edstomo`` module supports STEM/EDS tomographic reconstruction.  The code works as part of a larger workflow that can be seen in the examples in ncempy.data.  Please visit https://github.com/ercius/openNCEM/tree/master/ncempy/data/L2083-K-4-1 for a good start.

Contents
--------

Overview of contents with short description:

+---------------------------+--------------------------------------------------------------------+
| Module                    | Description                                                        |
+===========================+====================================================================+
| CharactersiticEmission    | X-ray fluorescence energies.                                       |
+---------------------------+--------------------------------------------------------------------+
| DoGenfire                 | Wrapper to run GENFIRE tomographic reconstructions.                |
+---------------------------+--------------------------------------------------------------------+
| Bruker                    | Move datasets from Bruker bcf files into emd files.                |
+---------------------------+--------------------------------------------------------------------+
| preprocess                | Tomography helper functions used before reconstruction.            |
+---------------------------+--------------------------------------------------------------------+
| postprocess               | Tomography helper functions used after reconstruction.             |
+---------------------------+--------------------------------------------------------------------+
