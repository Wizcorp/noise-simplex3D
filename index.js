var RandomNumberGenerator = require('math-random').RandomNumberGenerator;

/**
 * @classdesc 3-dimensional Simplex Noise
 * @class
 *
 * @author Brice Chevalier
 *
 * @param {object} params
 * @param {number} params.octaves
 * @param {number} params.amplitude
 * @param {number} params.frequency
 * @param {number} params.persistance
 * @param {number} params.base
 */

var grad = [
	[1, 1, 0],
	[-1, 1, 0],
	[1, -1, 0],
	[-1, -1, 0],
	[1, 0, 1],
	[-1, 0, 1],
	[1, 0, -1],
	[-1, 0, -1],
	[0, 1, 1],
	[0, -1, 1],
	[0, 1, -1],
	[0, -1, -1],
	[-1, 1, 1],
	[1, -1, 1],
	[1, 1, -1],
	[0, 0, 0]
];


function Simplex3D(params) {
	params = params || {};
	this.octaves = !params.octaves ? 1 : params.octaves;
	this.amplitude = !params.amplitude ? 1 : params.amplitude;
	this.frequency = !params.frequency ? 1 : params.frequency;
	this.persistance = !params.persistance ? 0.5 : Math.min(Math.max(params.persistance, 0), 1);

	// The scale is used to put the noise value in the interval [-amplitude / 2; amplitude / 2]
	this.scale = (this.persistance === 1) ? this.octaves * this.amplitude / 2 : (1 - this.persistance) / (1 - Math.pow(this.persistance, this.octaves)) * this.amplitude / 2;

	// The base is used to put the noise value in the interval [base; amplitude + base]
	this.base = (params.base || 0) + this.amplitude / 2;

	// initialize the permutation table
	this.seed(params.seed || 0);
}

/** Initialize noise's permutation table with provided seed
 *
 * @param {Number} seedNumber
 */
Simplex3D.prototype.seed = function (seedNumber) {
	var i;

	// reset permutation table
	var perm = this.perm = [];
	for (i = 0; i < 256; i++) perm.push(i);

	// randomly permute elements in table
	var random = new RandomNumberGenerator(seedNumber);
	for (i = 0; i < 256; i++) {
		var index = ~~(256 * random.next());
		// permute the two indexes
		var v = perm[i];
		perm[i] = perm[index];
		perm[index] = v;
	}

	// concat the table with itself to duplicate the permutations
	perm = perm.concat(perm);
};


var f3 = 1.0 / 3.0;
var g3 = 1.0 / 6.0;
Simplex3D.prototype.generateNoise = function (xin, yin, zin) {
	var perm = this.perm;

	var n0, n1, n2, n3; // Noise contributions from the four corners

	// Skew the input space to determine which simplex cell we're in
	var s = (xin + yin + zin) * f3; // Simple skew factor for 3D
	var i = Math.floor(xin + s);
	var j = Math.floor(yin + s);
	var k = Math.floor(zin + s);
	var t = (i + j + k) * g3;

	var x0 = i - t; // Unskew the cell origin back to (x,y,z) space
	var y0 = j - t;
	var z0 = k - t;
	x0 = xin - x0; // The x,y distances from the cell origin
	y0 = yin - y0;
	z0 = zin - z0;

	// For the 3D case, the simplex shape is a slightly irregular tetrahedron. // Determine which simplex we are in.
	var i1, j1, k1; // Offsets for second corner of simplex in (i,j,k) coords
	var i2, j2, k2; // Offsets for third corner of simplex in (i,j,k) coords
	if (x0 >= y0) {
		if (y0 >= z0) {
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0; // X Y Z order
		} else if (x0 >= z0) {
			i1 = 1;
			j1 = 0;
			k1 = 0;
			i2 = 1;
			j2 = 0;
			k2 = 1; // X Z Y order
		} else {
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 1;
			j2 = 0;
			k2 = 1; // Z X Y order
		}
	} else { // x0<y0
		if (y0 < z0) {
			i1 = 0;
			j1 = 0;
			k1 = 1;
			i2 = 0;
			j2 = 1;
			k2 = 1; // Z Y X order
		} else if (x0 < z0) {
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 0;
			j2 = 1;
			k2 = 1; // Y Z X order
		} else {
			i1 = 0;
			j1 = 1;
			k1 = 0;
			i2 = 1;
			j2 = 1;
			k2 = 0; // Y X Z order
		}
	}

	// A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z),
	// a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and // a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where // c = 1/6.
	var x1 = x0 - i1 + g3; // Offsets for second corner in (x,y,z) coords
	var y1 = y0 - j1 + g3;
	var z1 = z0 - k1 + g3;
	var x2 = x0 - i2 + 2.0 * g3; // Offsets for third corner in (x,y,z) coords
	var y2 = y0 - j2 + 2.0 * g3;
	var z2 = z0 - k2 + 2.0 * g3;
	var x3 = x0 - 1.0 + 3.0 * g3; // Offsets for last corner in (x,y,z) coords
	var y3 = y0 - 1.0 + 3.0 * g3;
	var z3 = z0 - 1.0 + 3.0 * g3;

	// Work out the hashed gradient indices of the three simplex corners
	var ii = i & 255;
	var jj = j & 255;
	var kk = k & 255;

	// Calculate the contribution from the three corners
	var t0 = 0.5 - x0 * x0 - y0 * y0 - z0 * z0;
	var t1 = 0.5 - x1 * x1 - y1 * y1 - z1 * z1;
	var t2 = 0.5 - x2 * x2 - y2 * y2 - z2 * z2;
	var t3 = 0.5 - x3 * x3 - y3 * y3 - z3 * z3;

	if (t0 < 0) {
		n0 = 0.0;
	} else {
		var gi0 = perm[ii + perm[jj + perm[kk]]] & 15;
		t0 *= t0;
		n0 = t0 * t0 * (grad[gi0][0] * x0 + grad[gi0][1] * y0 + grad[gi0][2] * z0);
	}

	if (t1 < 0) {
		n1 = 0.0;
	} else {
		var gi1 = perm[ii + i1 + perm[jj + j1 + perm[kk + k1]]] & 15;
		t1 *= t1;
		n1 = t1 * t1 * (grad[gi1][0] * x1 + grad[gi1][1] * y1 + grad[gi1][2] * z1);
	}

	if (t2 < 0) {
		n2 = 0.0;
	} else {
		var gi2 = perm[ii + i2 + perm[jj + j2 + perm[kk + k2]]] & 15;
		t2 *= t2;
		n2 = t2 * t2 * (grad[gi2][0] * x2 + grad[gi2][1] * y2 + grad[gi2][2] * z2);
	}

	if (t3 < 0) {
		n3 = 0.0;
	} else {
		var gi3 = perm[ii + 1 + perm[jj + 1 + perm[kk + 1]]] & 15;
		t3 *= t3;
		n3 = t3 * t3 * (grad[gi3][0] * x3 + grad[gi3][1] * y3 + grad[gi3][2] * z3);
	}

	// Add contributions from each corner to get the final noise value.
	// The result is scaled to return values in the interval [-1,1].
	return 32.0 * (n0 + n1 + n2 + n3);
};

// Complexity in O(o)
// with o the number of octaves
Simplex3D.prototype.getNoise = function (x, y, z) {
	var noise = 0;
	var amp = 1.0;

	for (var o = 0; o < this.octaves; o += 1) {
		noise += this.generateNoise(x, y, z) * amp;
		x *= this.frequency;
		y *= this.frequency;
		z *= this.frequency;
		amp *= this.persistance;
	}

	return noise * this.scale + this.base;
};

module.exports = Simplex3D;
