How does it work?
=================

Each piece of a phylogenetics projects is stored in a Handler. Multiple of the same
Handlers can be stored in a Handler subclass called HandlerContainer. The Project
class manages the complete mapping of the Handlers. Each Handler is completely decoupled
from all others. Two Handlers (and HandlerContainers) can be linked, which creates a
one-way soft reference between two Handlers. The main underpinning of a phylogenetics project
is the metadata attribute of Handlers. The API inits, updates, removes, and changes the
metadata underneath. When saving to disk, the metadata is written out, not the Python
API itself. This keeps the "dumping to disk" lightweight and static as the phylogenetics
API might rapidly change.
