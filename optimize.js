window._randNorm = null;
window.randomNormal = function () {
    // Box-Muller transform for normally distributed random numbers.
    // http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
    var f, u, v, s = 0.0;
    if (window._randNorm !== null &&
            typeof(window._randNorm) !== "undefined") {
        var tmp = window._randNorm;
        window._randNorm = null;
        return tmp;
    }
    while (s === 0.0 || s >= 1.0) {
        u = 2 * Math.random() - 1;
        v = 2 * Math.random() - 1;
        s = u * u + v * v;
    }
    f = Math.sqrt(-2 * Math.log(s) / s);
    window._randNorm = v * f;
    return u * f;
};

window.optimize = (function () {

    var optimize = {}, vector = {}, _q = function (x) {
        // The existential operator;
        return typeof(x) !== "undefined" && x !== null;
    };

    // ======= //
    //         //
    // VECTORS //
    //         //
    // ======= //

    vector.copy = function (x) {
        var y, i;
        if (typeof(x.length) === "undefined") return x;
        y = new Array();
        for (i = 0; i < x.length; i++) y[i] = x[i];
        return y;
    };

    vector.atleast_1d = function (x) {
        // Make sure that an object acts as an array (even if it's a scalar).
        if (typeof(x.length) === "undefined") {
            var tmp = new Array();
            tmp[0] = x;
            return tmp;
        }
        return x;
    };

    vector.range = function (a, b, c) {
        var xmin, xmax, dx, x, rng = new Array();
        if (typeof(b) === "undefined") {
            xmin = 0;
            xmax = a;
            dx = 1;
        } else if (typeof(c) === "undefined") {
            xmin = a;
            xmax = b;
            dx = 1;
        } else {
            xmin = a;
            xmax = b;
            dx = c;
        }
        for (x = xmin, i = 0; x < xmax; x += dx, i++) rng[i] = x;
        return rng;
    };

    vector.dot = function (a, b) {
        var i, result = 0.0;
        if (a.length !== b.length) throw "Size mismatch in vector.dot.";
        for (i = 0; i < a.length; i++) result += a[i] * b[i];
        return result;
    };

    vector.fmult = function (f, v) {
        var i, result = new Array();
        for (i = 0; i < v.length; i++) result[i] = f * v[i];
        return result;
    };

    vector.add = function (a, b) {
        var i, result = new Array();
        if (a.length !== b.length) throw "Size mismatch in vector.add.";
        for (i = 0; i < a.length; i++) result[i] = a[i] + b[i];
        return result;
    };

    vector.subtract = function (a, b) {
        var i, result = new Array();
        if (a.length !== b.length) throw "Size mismatch in vector.subtract.";
        for (i = 0; i < a.length; i++) result[i] = a[i] - b[i];
        return result;
    };

    vector.reduce = function (x) {
        var i, result = x[0];
        for (i = 1; i < x.length; i++)
            result = vector.add(result, x[i]);
        return result;
    };

    vector.take = function (x, ind) {
        // Re-order a vector.
        var i, result = new Array();
        for (i = 0; i < ind.length; i++) result[i] = x[ind[i]];
        return result;
    };

    vector.argsort = function (x) {
        // Argsort of an array using merge sort.
        if (typeof(x.length) === "undefined" || x.length === 1) return x;
        return vector._rec_argsort(vector.range(x.length), x);
    };

    vector._rec_argsort = function (inds, data) {
        // The recursive helper function for the argsort.
        var m, l, r;
        if (typeof(inds) === "undefined" || inds.length === 1) return inds;
        m = parseInt(inds.length / 2);
        l = inds.slice(0, m);
        r = inds.slice(m, inds.length);
        return vector._merge(vector._rec_argsort(l, data),
                             vector._rec_argsort(r, data), data);
    };

    vector._merge = function (l, r, data) {
        // Merging for use with a merge argsort.
        var result = new Array();
        while (l.length && r.length) {
            if (data[l[0]] <= data[r[0]]) result.push(l.shift());
            else result.push(r.shift());
        }
        while (l.length) result.push(l.shift());
        while (r.length) result.push(r.shift());
        return result;
    };

    // ============ //
    //              //
    // OPTIMIZATION //
    //              //
    // ============ //

    optimize._approx_fprime = function (x, f, ep) {
        // Approximate the N dimensional gradient of a scalar function using
        // forward finite difference.
        var i, f0 = f(x), grad = new Array();
        x = vector.atleast_1d(x);
        if (typeof(ep.length) === "undefined") {
            eps = [];
            for (i = 0; i < x.length; i++)
                eps.push(ep);
        } else if (ep.length === x.length) {
            eps = ep;
        } else throw "Size mismatch in _approx_fprime.";

        for (i = 0; i < x.length; i++) {
            x[i] += eps[i];
            grad[i] = (f(x) - f0) / eps[i];
            x[i] -= eps[i];
        }

        return grad;
    };

    optimize._approx_jacobian = function (x, f, ep) {
        // Approximate NxM dimensional gradient of the vector function
        // f using forward finite difference.
        var i, x0, f0 = $V(f(x)), grad = new Array();

        x = vector.atleast_1d(x);
        x0 = vector.copy(x);

        if (typeof(ep.length) === "undefined") {
            eps = [];
            for (i = 0; i < x.length; i++)
                eps.push(ep);
        } else if (ep.length === x.length) {
            eps = ep;
        } else throw "Size mismatch in _approx_fprime.";

        for (i = 0; i < x.length; i++) {
            x[i] = x0[i] + eps[i];
            grad[i] = $V(f(x)).subtract(f0).x(1.0 / eps[i]).elements;
            x[i] = x0[i];
        }

        return $M(grad).transpose();
    };

    optimize._max_abs_diff = function (x0, x) {
        var i, max = 0.0;
        for (i = 0; i < x.length; i++)
            max = Math.max(max, Math.abs(x0 - x[i]));
        return max;
    };

    optimize.fmin = function (func, x0, opts) {
        // Optimize a function using Nelder-Mead.
        var N, rho, chi, psi, sigma, sim, fsim, i, j, iterations;
        var nonzdelt, zdelt, x, fval;

        x0 = vector.atleast_1d(x0);
        N = x0.length;

        // Defaults.
        if (!_q(opts)) opts = {};
        xtol = _q(opts.xtol) ? opts.xtol : 1e-6;
        ftol = _q(opts.ftol) ? opts.ftol : 1e-6;
        maxiter = _q(opts.maxiter) ? opts.maxiter : 200 * N;

        // Magic numbers from `scipy`.
        rho = 1;
        chi = 2;
        psi = 0.5;
        sigma = 0.5;
        nonzdelt = 0.05;
        zdelt = 0.00025;

        sim = new Array();
        sim[0] = x0;
        fsim = new Array();
        fsim[0] = func(x0);

        for (i = 0; i < N; i++) {
            y = vector.copy(x0);
            if (y[i] !== 0.0) y[i] *= (1 + nonzdelt);
            else y[i] = zdelt;

            sim[i + 1] = y;
            fsim[i + 1] = func(y);
        }

        inds = vector.argsort(fsim);
        fsim = vector.take(fsim, inds);
        sim = vector.take(sim, inds);

        iterations = 0;

        // Constraint on function calls is needed.
        while (iterations < maxiter) {
            var xbar, xr, fxr, doshrink = false;
            iterations += 1
            // A break based on xtol needs to be included too.
            if (optimize._max_abs_diff(fsim[0], fsim.slice(1, fsim.length))
                    <= ftol)
                break;

            xbar = vector.fmult(1.0 / N,
                    vector.reduce(sim.slice(0, sim.length - 1)));
            xr = vector.add(vector.fmult(1 + rho, xbar),
                            vector.fmult(-rho, sim[sim.length - 1]));
            fxr = func(xr);

            if (fxr < fsim[0]) {
                var xe, fxe;
                xe = vector.add(vector.fmult(1 + rho * chi, xbar),
                            vector.fmult(-rho * chi, sim[sim.length - 1]));
                fxe = func(xe);
                if (fxe < fxr) {
                    sim[sim.length - 1] = xe;
                    fsim[fsim.length - 1] = fxe;
                } else {
                    sim[sim.length - 1] = xr;
                    fsim[fsim.length - 1] = fxr;
                }
            } else {
                if (fxr < fsim[fsim.length - 2]) {
                    sim[sim.length - 1] = xr;
                    fsim[fsim.length - 1] = fxr;
                } else {
                    var xc, fxc;
                    if (fxr < fsim[fsim.length - 1]) {
                        xc = vector.add(vector.fmult(1 + rho * psi, xbar),
                            vector.fmult(-rho * psi, sim[sim.length - 1]));
                        fxc = func(xc);
                        if (fxc < fsim[fsim.length - 1]) {
                            sim[sim.length - 1] = xc;
                            fsim[fsim.length - 1] = fxc;
                        } else {
                            doshrink = true;
                        }
                    } else {
                        xc = vector.add(vector.fmult(1 - psi, xbar),
                            vector.fmult(psi, sim[sim.length - 1]));
                        fxc = func(xc);
                        if (fxc < fsim[fsim.length - 1]) {
                            sim[sim.length - 1] = xc;
                            fsim[fsim.length - 1] = fxc;
                        } else {
                            doshrink = true;
                        }
                    }

                    if (doshrink) {
                        for (j = 1; j < N + 1; j++) {
                            sim[j] = vector.add(sim[0], vector.fmult(sigma,
                                            vector.subtract(sim[j], sim[0])));
                            fsim[j] = func(sim[j]);
                        }
                    }
                }
            }

            inds = vector.argsort(fsim);
            fsim = vector.take(fsim, inds);
            sim = vector.take(sim, inds);
        }

        x = sim[0];
        fval = fsim[0];

        if (iterations >= maxiter)
            console.log("Too many interations.", iterations);
        else
            console.log("Converged in", iterations, "iterations.");

        console.log("Function value =", fval);

        return x;
    };

    optimize.newton = function (fn, x0, opts) {
        // fn should return the chi vector (data - model) / sigma.
        // Also, this function uses _dumb-ass_ inversion. We suck.
        var N, ftol, maxiter, fprime, ep, alpha, J, JT, JTJ, diagJTJ,
            JTfx0, dx, fx, fx0, df;

        x0 = vector.atleast_1d(x0);
        chi0 = fn(x0);
        chi20 = vector.dot(chi0, chi0);

        // Defaults.
        if (!_q(opts)) opts = {};
        ftol = _q(opts.ftol) ? opts.ftol : 1e-10;
        ep = _q(opts.ep) ? opts.ep : 1.49e-8;  // Magic number from scipy.
        maxiter = _q(opts.maxiter) ? opts.maxiter : 200 * x0.length;
        fprime = _q(opts.fprime) ? opts.fprime : function (x) {
            return optimize._approx_fprime(x, fn, ep);
        };

        alpha = 1.0;

        for (i = 0; i < maxiter; i++) {
            J = optimize._approx_jacobian(x0, fn, ep);
            JT = J.transpose()
            JTJ = JT.x(J);
            diagJTJ = Matrix.Diagonal(JTJ.diagonal().elements);
            JTfx = JT.x($V(chi0));

            dx = JTJ.inv().x(JTfx);
            // dx = JTJ.add(diagJTJ.x(lambda)).inv().x(JTfx);

            x_best = vector.copy(x0);
            chi_best = vector.copy(chi0);
            chi2_best = chi20;

            for (n = 0; n <= 5; n++) {
                alpha = Math.pow(2, -n);

                x_try = vector.subtract(x0, dx.x(alpha).elements);
                chi_try = fn(x_try);
                chi2_try = vector.dot(chi_try, chi_try);

                if (chi2_try < chi2_best) {
                    x_best = x_try;
                    chi_best = chi_try;
                    chi2_best = chi2_try;
                }
            }

            dchi2 = chi20 - chi2_best;

            x0 = x_best;
            chi0 = chi_best;
            chi20 = chi2_best;

            if (dchi2 < 0.0) throw "Failure";

            if (i > 1 && dchi2 < ftol)
                break;
        }

        console.log("Converged after", i, "iterations.");

        return x0;
    };

    optimize.test = function () {
        var x, synth, data, chi, chi2, p0 = [10.5, 6.0], truth = [5.3, 3.0],
            p_newton, p_fmin;

        synth = function (x, noise) {
            var t, result = [], sky = x[0], sig = x[1], v, norm;
            v2 = sig * sig;
            norm = 1.0 / Math.sqrt(2 * Math.PI * v2);
            for (t = -10; t <= 10; t++)
                result.push(sky + Math.exp(-0.5 * t * t / v2) * norm +
                        noise * window.randomNormal());
            return result;
        };

        data = synth(truth, 0.1);

        chi = function (x) {
            return vector.subtract(data, synth(x, 0.0));
        };

        chi2 = function (x) {
            var f = chi(x);
            return vector.dot(f, f);
        };

        p_newton = optimize.newton(chi, p0);
        p_fmin = optimize.fmin(chi2, p0);

        console.log("truth:", truth);
        console.log("p_newton:", p_newton);
        console.log("p_fmin:", p_fmin);
    };

    optimize.vector = vector;

    return optimize;

})();
