optimize.js
===========

[This](http://en.wikipedia.org/wiki/Mathematical_optimization) kind of
optimization. Not [that](https://developers.google.com/closure/compiler/).

This was ported from the `scipy.optimize` and we seem to be about an order
of magnitude faster than the Python version but that hasn't been too well
tested.

Usage
-----

To optimize some dumb function, just run

```javascript
my_function = function (x) {
    return optimize.vector.dot(x, x);
};
xopt = optimize.fmin(my_function, [5.0, -3.4, 1.7, 16.3, 0.17]);
```

And this should say something like:

```
Converged in 349 iterations.
Function value = 3.954810202072493e-7
```

And `xopt` should end up being something like

```
[-0.0003993727670733724, -0.00027257793115893254, -0.0003811443958917447, 0.00012560283306980614, 0.00002523018782488204]
```

You can also include some options as follows:

```javascript
xopt = optimize.fmin(my_function, [5.0, -3.4, 1.7, 16.3, 0.17], {ftol: 1e-7, maxiter: 1000});
```

To Do
-----

* Write tests.
* Implement some better algorithms.
* Make some demos.