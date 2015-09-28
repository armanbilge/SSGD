/*
 * Integrator.java
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

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;

import java.util.HashMap;
import java.util.Map;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public abstract class Integrator extends AbstractModel {

    private final Map<ParameterValue,Double> probabilities = new HashMap<ParameterValue,Double>();

    public Integrator(String name) {
        super(name);
    }

    private static final class ParameterValue {
        private final int iState;
        private final double iTime;
        private final int jState;
        private final double jTime;
        private final double mu;

        private ParameterValue(int iState, double iTime, int jState, double jTime, double mu) {
            this.iState = iState;
            this.iTime = iTime;
            this.jState = jState;
            this.jTime = jTime;
            this.mu = mu;
        }

        @Override
        public boolean equals(final Object obj) {
            if (obj instanceof ParameterValue) {
                final ParameterValue other = (ParameterValue) obj;
                return ((iState == other.iState && iTime == other.iTime && jState == other.jState && jTime == other.jTime)
                        || (iState == other.jState && iTime == other.jTime && jState == other.iState && jTime == other.iTime))
                        && mu == other.mu;
            } else {
                return false;
            }
        }

        @Override
        public int hashCode() {
            int result;
            final long iTime = Double.doubleToLongBits(this.iTime);
            final long jTime = Double.doubleToLongBits(this.jTime);
            final int i = (int) (31 * iState + iTime ^ (iTime >>> 32));
            final int j = (int) (31 * jState + jTime ^ (jTime >>> 32));
            result = i ^ j;
            final long mu = Double.doubleToLongBits(this.mu);
            result = 31 * result + (int) (mu ^ (mu >>> 32));
            return result;
        }
    }


    public final double integratedProbability(int iState, double iTime, int jState, double jTime, double mu) {
        final ParameterValue value = new ParameterValue(iState, iTime, jState, jTime, mu);
        if (!probabilities.containsKey(value))
            probabilities.put(value, calculateIntegratedProbability(iState, iTime, jState, jTime, mu));
        return probabilities.get(value);
    }

    protected abstract double calculateIntegratedProbability(int iState, double iTime, int jState, double jTime, double mu);

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        fireModelChanged(o, i);
        probabilities.clear();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
        fireModelChanged(variable, i);
        probabilities.clear();
    }

    @Override
    protected void restoreState() {
        probabilities.clear();
    }

}
