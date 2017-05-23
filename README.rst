Crosstalk correction
--------------------
Python 3 implementation of `Donato et al.'s
2013 maximum impact estimation algorithm
<http://doi.org/10.1101/gr.153551.112>`_
for correcting crosstalk effects in pathway analysis.

Background
-----------

**Crosstalk**: Donato et al. (2013) use the term crosstalk to refer to the
effect that pathways exercise on each other (in pathway analysis methods
such as enrichment analysis, functional class scoring, and topology-based
methods) due to the presence of overlapping genes.

**Maximum impact estimation**: They developed a correction method called
maximum impact estimation that takes into account overlaps between pathways.
The approach infers an underlying pathway impact matrix where each gene
only contributes to one pathway using an expectation maximization technique.

**PathCORE**: The crosstalk correction method is used in the PathCORE software,
a hypothesis generation tool that identifies co-occurring pathways from the
results of an unsupervised analysis of transcriptomic data. Due to confusion
around the term "crosstalk," we refer to this procedure as "overlap-correction"
in the PathCORE software and paper.

Installation
----------------
Please use Python 3.4 or higher.
To install the current PyPI version (recommended), run::

    pip install crosstalk-correction

For the latest GitHub version, run::

    pip install git+https://github.com/kathyxchen/crosstalk-correction.git#egg=crosstalk-correction

Package contents
----------------

=======================
crosstalk_correction.py
=======================
crosstalk_correction.py contains the implementation of the crosstalk
correction procedure. The ``crosstalk_correction`` function wraps
the maximum impact estimation algorithm (in ``maximum_impact_estimation``)
and reduces the number of pre/post-processing steps required to
run/interpret the results of ``maximum_impact_estimation``.

We recommend that ``crosstalk_correction`` be used directly, but
the functions called within ``crosstalk_correction`` have also
been documented so that users can determine how to best address
their use case.

Acknowledgements
----------------
This work was supported by the Penn Institute for Bioinformatics
