/*
 * NegativeMultivariateFunction.java
 *
 * SSGD: Serially-Sampled Genome Demographics
 *
 * Copyright (c) 2015 Arman Bilge <armanbilge@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.compevol.ssgd;

import dr.math.MultivariateFunction;

/**
 * @author Arman Bilge
 */
public class NegativeMultivariateFunction implements MultivariateFunction {

    private final MultivariateFunction f;

    public NegativeMultivariateFunction(final MultivariateFunction f) {
        this.f = f;
    }

    @Override
    public double evaluate(double[] doubles) {
        return - f.evaluate(doubles);
    }

    @Override
    public int getNumArguments() {
        return f.getNumArguments();
    }

    @Override
    public double getLowerBound(int i) {
        return f.getLowerBound(i);
    }

    @Override
    public double getUpperBound(int i) {
        return f.getUpperBound(i);
    }

}
