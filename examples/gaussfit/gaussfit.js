window.model = (function () {

    var Model = function (a, noise) {
        var i;
        this.amp = a[0];
        this.mu = a[1];
        this.sig2 = a[2] * a[2];

        this.order = a.length - 3;
        this.a = [];
        for (i = 0; i < this.order; i++) this.a.push(a[i + 3]);

        if (typeof(noise) === "undefined") this.noise = 0.0;
        else this.noise = noise;
    };

    Model.prototype.predict = function (x, noise) {
        // Calculate the model prediction at a set of x points.
        var norm, i, j, result = [];
        if (typeof(noise) === "undefined") noise = 0.0;
        else noise = this.noise;
        x = optimize.vector.atleast_1d(x);
        norm = this.amp / Math.sqrt(2 * Math.PI * this.sig2);

        // Initialize the prediction with the Gaussian.
        for (i = 0; i < x.length; i++) {
            var diff = x[i] - this.mu;
            result.push(norm * Math.exp(-0.5 * diff * diff / this.sig2));
        }

        // Add in the polynomial background and optional noise.
        for (j = 0; j < this.order; j++) {
            for (i = 0; i < x.length; i++) {
                result[i] += this.a[j] * Math.pow(x[i], j)
                             + noise * randomNormal();
            }
        }

        return result;
    };

    Model.prototype.line = function (xaxis) {
        var i, data = [], yaxis = this.predict(xaxis);
        for (i = 0; i < yaxis.length; i++) {
            data.push({x: xaxis[i], y: yaxis[i]});
        }
        return data;
    };

    var fake = function (order) {
        var i, truth = [-2 * Math.random(), 0.1 + 0.01 * Math.random(),
            0.2 * (2 * Math.random() - 1)];
        for (i = 0; i < order; i++) truth.push(2 * randomNormal());
        return new Model(truth);
    };

    return {Model: Model,
        fake: fake};

})();

window.Dataset = (function () {

    var Dataset = function (ndata, noise) {
        var x;
        this.ndata = ndata;
        this.noise = noise;

        // Generate the x-axis for the line plots.
        this.xrange = [-1, 1];
        this.xaxis = []
        for (x = xrange[0]; x <= xrange[1];
                                    x += (xrange[1] - xrange[0]) / 500.0) {
            this.xaxis.push(x);
        }

        this.gen_data();
    };

    Dataset.prototype.gen_data = function () {
        var i;
        // Generate the fake data.
        this.truth = model.fake(6);
        this.x = [];
        this.yerr = [];
        for (i = 0; i < this.ndata; i++) {
            this.x.push((this.xrange[1] - this.xrange[0]) * Math.random()
                + this.xrange[0]);
            this.yerr.push(this.noise);
        }
        this.y = this.truth.predict(this.x, this.noise);

        this.data = [];
        for (i = 0; i < this.x.length; i++) {
            this.data.push({x: this.x[i], y: this.y[i], yerr: this.yerr[i]});
        }

        return this.data;
    };

    Dataset.prototype.line = function () {
        var i, data = [], yaxis = this.truth.predict(this.xaxis);
        for (i = 0; i < yaxis.length; i++) {
            data.push({x: this.xaxis[i], y: yaxis[i]});
        }
        return data;
    };

    Dataset.prototype.chi = function (p) {
        var i, m = new Model(p), ym, chi;
        ym = m.predict(this.x);
        for (i = 0; i < ym.length; i++) {
            chi.push((this.y[i] - ym[i]) / this.yerr[i]);
        }
        return chi;
    };

    Dataset.prototype.chi2 = function (p) {
        var i, chi2 = 0.0, chi = this.chi(p);
        for (i = 0; i < chi.length; i++) {
            chi2 += chi[i] * chi[i];
        }
        return chi2;
    };

    return Dataset;

})();
