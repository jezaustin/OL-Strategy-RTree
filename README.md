OL-Strategy-RTree
=================

A strategy to maintain a spatial index on a vector layer.

Operating on a feature's neighbours is impractical for a large vector layer.
A spatial index made available client side, within OpenLayers, would help.

At the moment, this is implemented as a strategy which users initialize with
event handlers to be called within a context from which the index is available
(this.Tree). The index is maintained as features are added and removed.
