Examples
========

Assembly
--------

Simulating DNA assembly between SeqRecords without a protocol.

Gibson Assembly with primer generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: gibson_assembly.py

Cloning via restriction digest and ligation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: clone_assembly.py

Protocols
---------

Executing a DNA assembly while accumulating a `Protocol` object for writing
protocol instructions, plate layouts, DNA files (`gb`, `fa`) and protocol inputs.

Gibson Assembly of a plasmid library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: gibson_library.py

Gibson Assembly of a single plasmid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: gibson_plasmid.py

Golden Gate Assembly of a single plasmid
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: goldengate_plasmid.py

Combinatorial Golden Gate Assembly: all valid plasmid combinations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: goldengate_combinatorial.py

Combinatorial Golden Gate Assembly with bins: all plasmid combinations between bins
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. literalinclude:: goldengate_combinatorial_bins.py
