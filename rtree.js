/*jslint continue: true, plusplus: true, vars: true, white: true, browser: true */
var OpenLayers;

var RStarTree = function RTree(M, m, p) {
    // beckman, kriegel, schneider, seeger 1990
    //http://dbs.mathematik.uni-marburg.de/publications/myPapers/1990/BKSS90.pdf
    'use strict';
    var Root, Branch, Leaf,

    // abstract node type, base of Leaf & Branch
    Node = function NodeFactory(parent) {
	var self = [],

	sortAlongAxis = function(axis) {
	    return function sortAlongAxis(a, b) {
		return (a.bounds[axis[0]] === b.bounds[axis[0]]) ?
		    a.bounds[axis[1]] - b.bounds[axis[1]] :
		    a.bounds[axis[0]] - b.bounds[axis[0]];
	    };
	},
	splitAxisScore = function splitAxisScore(dims) {
	    self.sort(sortAlongAxis(dims.split(0, 2)));

	    var S = 0, k = m,
	    getDim = function (dim) {
		return function getDim(x) {
		    return x.bounds[dim];
		};
	    },
	    A = self.slice(0, m),
	    B = self.slice(m),
	    minA1 = self[0].bounds[dims[0]], // due to sort
	    maxA1 = Math.max.apply(Math, A.map(getDim(dims[1]))),
	    minA2 = Math.min.apply(Math, A.map(getDim(dims[2]))),
	    maxA2 = Math.max.apply(Math, A.map(getDim(dims[3]))),
	    minB1 = B[0].bounds[dims[0]],
	    maxB1 = Math.max.apply(Math, B.map(getDim(dims[1]))),
	    minB2 = Math.min.apply(Math, B.map(getDim(dims[2]))),
	    maxB2 = Math.max.apply(Math, B.map(getDim(dims[3])));

	    while (true) { //k <= self.length - m) {
		// add boundary(A) + boundary(B)
		// (divided by 2)
		S += maxA1 - minA1 + maxA2 - minA2;
		S += maxB1 - minB1 + maxB2 - minB2;

		k++;
		if (k > self.length - m) {
		    break;
		}
		//A = self.slice(0, k); // A is not referred to!
		B = self.slice(k);
		//minA1 unchanged due to sort
		if (maxA1 < self[k-1].bounds[dims[1]]) {
		    maxA1 = self[k-1].bounds[dims[1]];
		}
		if (minA2 > self[k-1].bounds[dims[2]]) {
		    minA2 = self[k-1].bounds[dims[2]];
		}
		if (maxA2 < self[k-1].bounds[dims[3]]) {
		    maxA2 = self[k-1].bounds[dims[3]];
		}
		minB1 = B[0].bounds[dims[0]]; // due to sort
		if (maxB1 === self[k-1].bounds[dims[1]]) {
		    maxB1 = Math.max.apply(Math, B.map(getDim(dims[1])));
		}
		if (minB2 === self[k-1].bounds[dims[2]]) {
		    minB2 = Math.min.apply(Math, B.map(getDim(dims[2])));
		}
		if (maxB2 === self[k-1].bounds[dims[3]]) {
		    maxB2 = Math.max.apply(Math, B.map(getDim(dims[3])));
		}
	    }
	    return S;
	},

	distribute = function distribute(targets, axis) {
	    self.sort(sortAlongAxis(axis.major));
	    // minimise overlap
	    var k, i,
	    getBounds = ('Branch' === self.type) ?
		function getBoundsNode(node) {
		    return node.bounds;
		} : function getBoundsFeature(feature) {
		    return feature.geometry.getBounds();
		},
	    findMins = function findMins(test, subset) {
		var val, minVal = Infinity, min_ks,
		k = subset ? subset.shift() : m,

		A = self.slice(0, k),
		B = self.slice(k),
		boundA = {}, boundB = {},
		getDim = function(dim) {
		    return function getDim(node) {
			return getBounds(node)[dim];
		    };
		};

		boundA[axis.major[0]] = A[0].bounds[axis.major[0]];
		boundA[axis.major[1]] = Math.max.apply(Math, A.map(getDim(axis.major[1])));
		boundA[axis.minor[0]] =	Math.min.apply(Math, A.map(getDim(axis.minor[0])));
		boundA[axis.minor[1]] =	Math.max.apply(Math, A.map(getDim(axis.minor[1])));

		boundB[axis.major[0]] = B[0].bounds[axis.major[0]];
		boundB[axis.major[1]] = Math.max.apply(Math, B.map(getDim(axis.major[1])));
		boundB[axis.minor[0]] =	Math.min.apply(Math, B.map(getDim(axis.minor[0])));
		boundB[axis.minor[1]] =	Math.max.apply(Math, B.map(getDim(axis.minor[1])));

		while (true) {
		    val = test(boundA, boundB);
		    if (val < minVal) {
			min_ks = [k];
		    } else if (val === minVal) {
			min_ks.push(k);
		    }


		    k = subset ? subset.shift() : k+1;
		    if (undefined === k || k > self.length - m) {
			break;
		    }

		    // update bound[A|B]
		    if (boundA[axis.major[1]] < self[k-1].bounds[axis.major[1]]) {
			boundA[axis.major[1]] = self[k-1].bounds[axis.major[1]];
		    }
		    if (boundA[axis.minor[0]] > self[k-1].bounds[axis.minor[0]]) {
			boundA[axis.minor[0]] = self[k-1].bounds[axis.minor[0]];
		    }
		    if (boundA[axis.minor[1]] < self[k-1].bounds[axis.minor[1]]) {
			boundA[axis.minor[1]] = self[k-1].bounds[axis.minor[1]];
		    }

		    B = self.slice(k);
		    boundB[axis.major[0]] = B[0].bounds[axis.major[0]];
		    boundB[axis.major[1]] =
			Math.max.apply(Math, B.map(getDim(axis.major[1])));
		    boundB[axis.minor[0]] =
			Math.min.apply(Math, B.map(getDim(axis.minor[0])));
		    boundB[axis.minor[1]] =
			Math.max.apply(Math, B.map(getDim(axis.minor[1])));
		}
		return min_ks;
	    },
	    mins = findMins(function overlap(A, B) {
		return Math.max(0, Math.min(A.top, B.top) - Math.max(A.bottom, B.bottom)) *
		    Math.max(0, Math.min(A.right, B.right) - Math.min(A.left, B.left));
	    });
	    if (1 < mins.length) {
		mins = findMins(function area(A, B) {
		    return (A.top - A.bottom) * (A.right - A.left) +
			(B.top - B.bottom) * (B.right - B.left);
		}, mins);
	    }

	    k = mins.shift();

	    //distribute 0 ... k-1 in targets[0]; k ... in targets[1]
	    i = k;
	    while (i--) {
		targets[0].addElement(self[i]);
		targets[0].bounds.extend(getBounds(self[i]));
	    }
	    i = self.length;
	    while (k < i--) {
		targets[1].addElement(self[i]);
		targets[1].bounds.extend(getBounds(self[i]));
	    }
	},

	split = function split() {
	    // replace self's elements with two nodes;
	    // distribute self's elements between those two nodes.
	    var targets = ('Branch' === self.type) ?
		[new Branch(parent), new Branch(parent)] :
		[new Leaf(parent), new Leaf(parent)],
	    //replacement = new Branch(parent),

	    newRoot = (null === parent),

	    // choose split axis
	    splitAxis = (splitAxisScore(['left', 'right', 'bottom', 'top']) <
			 splitAxisScore(['bottom', 'top', 'left', 'right'])) ?
		{
		    major: ['left', 'right'],
		    minor: ['bottom', 'top']
		} : {
		    major: ['bottom', 'top'],
		    minor: ['left', 'right']
		};
	    // choose distribution along axis
	    distribute(targets, splitAxis);

	    if (newRoot) {
		parent = new Branch(null);
	    } else {
		parent.removeBranch(self);
	    }
	    parent.addBranch(targets[0]);
	    parent.addBranch(targets[1]);
	    if (newRoot) {
		Root = parent;
	    }
	},

	reInsert = function reInsert() {
	    // sort self by centroid distance
	    var boundsCentroid = self.bounds.getCentroid(), centroid,
	    i = self.length, reInsertions = [], reInsertion,
	    getCentroid = ('Branch' === self.type) ?
		function getCentroidNode(node) {
		    return node.bounds.getCentroid();
		} : function getCentroidFeature(feature) {
		    return feature.geometry.getCentroid();
		};
	    
	    while (i--) {
		centroid = getCentroid(self[i]);
		self[i].distCentroid2 =
		    Math.pow(centroid.x - boundsCentroid.x, 2) +
		    Math.pow(centroid.y - boundsCentroid.y, 2);
	    }
	    self.sort(function sortByDistance(A, B) {
		return B.distCentroid2 - A.distCentroid2;
	    });
	    
	    // remove p (0.3*M) most distant entries
	    i = p;
	    while (i--) {
		reInsertions.push(self.shift());
	    }
	    // refresh self.bounds

	    // reinsert them in order, closest first
	    reInsertion = reInsertions.pop();
	    while (reInsertion) {
		Root.insert(reInsertion);
		// Root.insert iterates down to Leaf level.
		// TODO: do a non-recursive parent.insert(reInsertion) instead?
		// OR: reinsert every feature in reInsertion separately?
		reInsertion = reInsertions.pop();
	    }
	};

	self.addElement = function addElement(element) {
	    self.push(element);
	    if (M < self.length) {
		if (self.reInsertGuard) {
		    split();
		} else {
		    self.reInsertGuard = true;
		    reInsert();
		    self.reInsertGuard = false;
		}
	    }
	};

	self.bounds = new OpenLayers.Bounds();

	return self;
    };
    Branch = function BranchFactory(parent) {
	var self = new Node(parent);

	self.insert = function ChooseBranch(feature) {
	    var isFeature = !feature.hasOwnProperty('bounds'),
	    getBounds = isFeature ?
		function getBoundsFeature(feature) {
		    return feature.geometry.getBounds();
		} : function getBoundsNode(node) {
		    return node.bounds.clone();
		};
    
	    // cost functions:
	    function area(i) {
		// return i's total area, including feature
		// (self[i] could be leaf or branch)
		var expanded = getBounds(feature);
		expanded.extend(self[i].bounds);
		return expanded.getWidth() * expanded.getHeight();
	    }
	    function enlargement(i) {
		// return enlargement i'th rectangle would need to accommodate feature
		// (self[i] could be leaf or branch)
		var expanded = getBounds(feature);
		expanded.extend(self[i].bounds);
		return (expanded.getHeight() - self.getHeight()) * self.getWidth() +
		    (expanded.getWidth() - self.getWidth()) * expanded.getHeight();
	    }
	    function overlap(i) {
		// return Sum A(E_i n E_j)
		// (self[i] are all leaves)
		var j = self.length, sum = 0, intersect;
		while (j--) {
		    if (i === j) {
			continue;
		    }
		    if (!self[i].bounds.intersectsBounds(self[j].bounds)) {
			continue;
		    }
		    intersect = new OpenLayers.Bounds();
		    intersect.extend(new OpenLayers.LonLat(
			Math.max(self[i].bounds.left, self[j].bounds.left),
			Math.max(self[i].bounds.bottom, self[j].bounds.bottom)
		    ));
		    intersect.extend(new OpenLayers.LonLat(
			Math.min(self[i].bounds.right, self[j].bounds.right),
			Math.min(self[i].bounds.top, self[j].bounds.top)
		    ));
		    sum += intersect.getWidth() * intersect.getHeight();
		}
		return sum;
	    }

	    var findMins = function findMins(test, subset) {
		var val, minVal = Infinity, minNodes,
		i = subset ? self.length : subset.length;
		while (i--) {
		    val = test(subset ? subset[i] : i);
		    if (val < minVal) {
			minVal = val;
			minNodes = [subset ? subset[i] : i];
		    } else if (val === minVal) {
			minNodes.push(subset ? subset[i] : i);
		    }
		}
		return minNodes;
	    },

	    mins;
	    // assume self non-empty!
	    if ('Branch' === self[0].type) {
		// minimum area cost
		mins = findMins(enlargement, null);
		if (1 < mins.length) {
		    mins = findMins(area, mins);
		}

	    } else {
		if (!isFeature) {
		    self.addBranch(feature);
		    self.bounds.extend(feature.bounds);
		    return;
		}

		// minimum overlap cost
		mins = findMins(overlap, null);
		if (1 < mins.length) {
		    mins = findMins(enlargement, mins);
		    if (1 < mins.length) {
			mins = findMins(area, mins);
		    }
		}
	    }

	    self[mins[0]].insert(feature);
	    self.bounds.extend(feature.geometry.getBounds());
	};

	self.search = function search(geometry) {
	    var results = [], i = self.length;
	    while (i--) {
		if (geometry.intersects(self[i].bounds)) {
		    Array.push.apply(results, self[i].search(geometry));
		}
	    }
	};

	self.addBranch = self.addElement;
	self.type = 'Branch';
	return self;
    };
    Leaf = function LeafFactory(parent) {
	var self = new Node(parent);

	self.insert = function insert(feature) {
	    self.addElement(feature);
	    self.bounds.extend(feature.geometry.getBounds());
	};

	self.search = function search(geometry) {
	    return self.filter(function featureIntersects(feature) {
		return geometry.intersects(feature.geometry);
	    });
	};

	self.type = 'Leaf';
	return self;
    };

    Root = new Leaf();

    return {
	search: function(geometry) {
	    Root.search(geometry);
	},
	insert: function(feature) {
	    Root.insert(feature);
	}
    };
};
// Production steps of ECMA-262, Edition 5, 15.4.4.19
// Reference: http://es5.github.com/#x15.4.4.19
if (!Array.prototype.map) {
    Array.prototype.map = function(callback, thisArg) {
	'use strict';

	var T, A, k;
 
	if (this == null) {
	    throw new TypeError(" this is null or not defined");
	}
 
	// 1. Let O be the result of calling ToObject passing the |this| value as the argument.
	var O = Object(this);
 
	// 2. Let lenValue be the result of calling the Get internal method of O with the argument "length".
	// 3. Let len be ToUint32(lenValue).
	var len = O.length >>> 0;
 
	// 4. If IsCallable(callback) is false, throw a TypeError exception.
	// See: http://es5.github.com/#x9.11
	if (typeof callback !== "function") {
	    throw new TypeError(callback + " is not a function");
	}
 
	// 5. If thisArg was supplied, let T be thisArg; else let T be undefined.
	if (thisArg) {
	    T = thisArg;
	}
 
	// 6. Let A be a new array created as if by the expression new Array(len) where Array is
	// the standard built-in constructor with that name and len is the value of len.
	A = new Array(len);
 
	// 7. Let k be 0
	k = 0;
 
	// 8. Repeat, while k < len
	while(k < len) {
 
	    var kValue, mappedValue;
 
	    // a. Let Pk be ToString(k).
	    //   This is implicit for LHS operands of the in operator
	    // b. Let kPresent be the result of calling the HasProperty internal method of O with argument Pk.
	    //   This step can be combined with c
	    // c. If kPresent is true, then
	    if (k in O) {
 
		// i. Let kValue be the result of calling the Get internal method of O with argument Pk.
		kValue = O[ k ];
 
		// ii. Let mappedValue be the result of calling the Call internal method of callback
		// with T as the this value and argument list containing kValue, k, and O.
		mappedValue = callback.call(T, kValue, k, O);
 
		// iii. Call the DefineOwnProperty internal method of A with arguments
		// Pk, Property Descriptor {Value: mappedValue, : true, Enumerable: true, Configurable: true},
		// and false.
 
		// In browsers that support Object.defineProperty, use the following:
		// Object.defineProperty(A, Pk, { value: mappedValue, writable: true, enumerable: true, configurable: true });
 
		// For best browser support, use the following:
		A[ k ] = mappedValue;
	    }
	    // d. Increase k by 1.
	    k++;
	}
 
	// 9. return A
	return A;
    };
}

if (!Array.prototype.filter) {
    Array.prototype.filter = function(fun /*, thisp */) {
	"use strict";
 
	if (this == null)
	    throw new TypeError();
 
	var t = Object(this);
	var len = t.length >>> 0;
	if (typeof fun != "function")
	    throw new TypeError();
 
	var res = [];
	var thisp = arguments[1];
	for (var i = 0; i < len; i++) {
	    if (i in t) {
		var val = t[i]; // in case fun mutates this
		if (fun.call(thisp, val, i, t))
		    res.push(val);
	    }
	}
 
	return res;
    };
}


var protoRTree = function RTree(wSplit) {
    'use strict';
    var //wSplit = 7,
    //wMerge = 3, // would use if we supported deleting features!
    Branch = function BranchFactory(parent) {
	return {
	    parent: parent,
	    branches: [],
	    leaves: [],
	    bounds: new OpenLayers.Bounds(),
	    insert: function insert(feature) {
		if (this.branches.length + this.leaves.length < wSplit) {
		    this.leaves.push(feature);
		    this.bounds.extend(feature.geometry.getBounds());
		    return;
		}

		this.branches[this.chooseBranch(feature)].insert(feature);
		this.bounds.extend(feature.geometry.getBounds());
	    },
	    boundsGrowth: function boundsGrowth(feature) {
		var expanded = feature.geometry.getBounds();
		expanded.extend(this.bounds);
		return this.bounds.getWidth() * (expanded.getHeight() - this.bounds.getHeight()) +
		    expanded.getHeight() * (expanded.getWidth() - this.bounds.getWidth());
	    },
	    chooseBranch: function chooseBranch(feature) {
		// loop over branches & leaves, choose most suitable
		// if choose a leaf, handle insert here!
		var expansion, minExpansion = Infinity, minBranch, minLeaf = null,
		i = this.branches.length;
		while (i--) {
		    expansion = this.branches[i].boundsGrowth(feature);
		    if (minExpansion > expansion) {
			minExpansion = expansion;
			minBranch = i;
		    }
		}

		// check leaves
		i = this.leaves.length;
		var bound;
		while (i--) {
		    // compare size of expansion to minExpansion
		    bound = feature.geometry.getBounds();
		    bound.extend(this.leaves[i].geometry);
		    expansion = bound.getHeight() * bound.getWidth();
		    if (minExpansion > expansion) {
			minExpansion = expansion;
			minLeaf = i;
		    }
		}
		if (null === minLeaf) {
		    return i;
		}

		// replace leaves[minLeaf] into a branch with 2 leaves!
		var branch = new Branch(this);
		branch.insert(this.leaves.splice(minLeaf, 1).pop());
		this.branches.pop(branch);
		return this.branches.length - 1;
	    },

	    findLeavesWithin: function findLeavesWithin(geom) {
		if (!this.bounds.toGeometry().intersects(geom)) {
		    return [];
		}
		var features = [],
		i = this.branches.length;
		while (i--) {
		    Array.push.apply(features, this.branches[i].findLeavesWithin(geom));
		}
		i = this.leaves.length;
		while (i--) {
		    if (this.leaves[i].geometry.intersects(geom)) {
			features.push(this.leaves[i]);
		    }
		}
		return features;
	    },
	    destroy: function destroy() {
		var i = this.branches.length;
		while (i--) {
		    this.branches[i].destroy();
		    delete this.branches[i];
		}
		// no need to actually destroy leaves (features)
		delete this.bounds;
	    }
	};
    },

    /*
    // Example of how to use findNeighboursWithin:
    findNeighbours = function findNeighbours(feature, radius) {
	// find points within a radius of feature
	var circle = OpenLayers.Geometry.Polygon.createRegularPolygon(feature.geometry.getCentroid(), radius, 20);
	return root.findNeighboursWithin(circle);
    }
    */
    root = new Branch(null);
    return root;
};

/*
usage:
new OpenLayers.Layer.Vector({
    protocol: blah...,
    // etc...
    strategies: [
	new OpenLayers.Strategy.RTree({
            autoActivate: false,
	    indexedListeners: [
		{
		    source: hoverControl,
		    event: 'featureunhighlighted',
		    listener: onFeatureUnHighlighted
		},{
		    source: map,
		    event: 'moveEnd',
		    listener: triggerUpdateStyleByDensity
		}//,etc.
	    ]
	})
    ]
});
// initialise hoverControl etc., then .activate() the strategy
*/



OpenLayers.Strategy.RTree = OpenLayers.Class(OpenLayers.Strategy, {
    RTree: null,
    activate: function activate() {
	'use strict';
        var activated = OpenLayers.Strategy.prototype.activate.call(this);
	if (activated) {
	    this.clearIndex();
	    this.createIndex({features: this.layer.features});

	    // events to keep up the index
	    this.layer.events.on({
		"beforefeaturesadded": this.createIndex,
		"featuresremoved": this.clearIndex,
		scope: this
	    });

	    // events to consume the index
	    var i = this.indexedListeners.length;
	    while (i--) {
		this.indexedListeners[i].source.events.register(
		    this.indexedListeners[i].event, this,
		    this.indexedListeners[i].listener);
	    }	     
	}
	return activated;
    },
    deactivate: function deactivate() {
	'use strict';
	var deactivated = OpenLayers.Strategy.prototype.deactivate.call(this);
	if (deactivated) {
	    this.clearIndex();
	    this.layer.events.un({
		"beforefeaturesadded": this.createIndex,
		"featuresremoved": this.clearIndex,
		scope: this
	    });

	    var i = this.indexedListeners.length;
	    while (i--) {
		this.indexedListeners[i].source.events.unregister(
		    this.indexedListeners[i].event, this,
		    this.indexedListeners[i].listener);
	    }
	}
	return deactivated;
    },
    createIndex: function createIndex(event) {
	'use strict';
	var i = event.features.length;
	while (i--) {
	    this.RTree.insert(event.features[i]);
	}
    },
    clearIndex: function clearIndex() {
	'use strict';
	if (this.RTree) {
	    this.RTree.destroy();
	}
	delete this.RTree;
	this.RTree = new RTree(7);
    },
    CLASS_NAME: "OpenLayers.Strategy.RTree"
});
