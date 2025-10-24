ncempy.eval
===========

The ``ncempy.eval`` module contains evaluation routines to address various data evaluation or simulation problems. They are made combining the algorithms found in ``ncempy.algo`` acting as bulding blocks for this higher layer of abstraction. An evaluation contained in ``ncempy.eval`` targets a specific problem and should not be reused/imported into another evaluation.

Contents
--------

Overview of contents with short description:

+--------------------+--------------------------------------------------------------------+
| Module             | Description                                                        |
+====================+====================================================================+
| ring_diff          | Evaluate ring diffraction patterns.                                |
+====================+====================================================================+
| line_profile       | Calculate line profiles along arbitrary angles in 2D images.       |
+====================+====================================================================+
| multicorr          | Advanced correlation of images to calculate shifts.                |
+====================+====================================================================+
| stack_align        | Alignment of stack of images                                       |
+--------------------+--------------------------------------------------------------------+

