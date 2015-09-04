/*
 * HKYSkylineIntegrator.java
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

import dr.evolution.coalescent.PiecewiseConstantPopulation;
import dr.evomodel.coalescent.PiecewisePopulationModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman Bilge
 */
public class HKYSkylineIntegrator extends Integrator {

    private final HKY hky;
    private final FrequencyModel frequencyModel;
    private final PiecewisePopulationModel populationModel;

    private boolean betaKnown = false;
    private double beta;
    private double kappa;

    public HKYSkylineIntegrator(final HKY hky, final PiecewisePopulationModel populationModel) {
        super("HKYSkylineIntegrator");
        this.hky = hky;
        addModel(hky);
        frequencyModel = hky.getFrequencyModel();
        this.populationModel = populationModel;
        addModel(populationModel);
    }

    private void calculateBeta() {
        kappa = hky.getKappa();
        final double freqA = frequencyModel.getFrequency(0);
        final double freqC = frequencyModel.getFrequency(1);
        final double freqG = frequencyModel.getFrequency(2);
        final double freqT = frequencyModel.getFrequency(3);
        final double freqR = freqA + freqG;
        final double freqY = freqC + freqT;
        beta = 1.0 / (2 * (freqR * freqY + kappa * (freqA * freqG + freqC * freqT)));
        betaKnown = true;
    }

    @Override
    protected double calculateIntegratedProbability(final int iState, final double iTime, final int jState, final double jTime, final double mu) {

        if (!betaKnown)
            calculateBeta();

        final double tau = Math.abs(iTime - jTime);

        final H H;
        if (iState % 2 == jState % 2) { // transition
            H = new H() {
                @Override
                public double apply(final double t, final double N) {
                    return transitionH(t, N, iState, jState, tau, mu);
                }
            };
        } else { // transversion
            H = new H() {
                @Override
                public double apply(final double t, final double N) {
                    return transversionH(t, N, jState, tau, mu);
                }
            };
        }

        return integrateIntervals(H, Math.max(iTime, jTime));
    }

    private double integrateIntervals(final H H, final double start) {

        final PiecewiseConstantPopulation df = (PiecewiseConstantPopulation) populationModel.getDemographicFunction();
        final int m = df.getNumArguments();

        int k;
        double current;
        for (k = 0, current = 0; current <= start; current += df.getEpochDuration(k++));
        double previous = start;

        double g = 1.0;
        double integratedP = 0;
        for (int i = k; i < m; ++i) {

            final double N = df.getEpochDemographic(i - 1);
            integratedP += g * Math.exp(previous / N) * (H.apply(current, N) - H.apply(previous, N));

            g *= Math.exp(-(current - previous) / N);

            previous = current;
            current += df.getEpochDuration(i);

        }

        final double N = df.getEpochDemographic(m - 1);
        integratedP -= g * Math.exp(previous / N) * H.apply(previous, N);

        if (Double.isNaN(integratedP))
            return 0.0;

        return integratedP;

    }

    private interface H {
        double apply(double t, double N);
    }

    private double transitionH(final double t, final double N, final int i, final int j, final double tau, final double mu) {

        final int ihat = (i + 2) % 4;
        final int pm = i == j ? -1 : 1;

        final double betamu = beta * mu;
        final double twobetamuN = 2 * betamu * N;
        final double mbetamutwotptau = -betamu * (2*t + tau);

        final double freqi = frequencyModel.getFrequency(i);
        final double freqj = frequencyModel.getFrequency(j);
        final double freqihat = frequencyModel.getFrequency(ihat);
        final double freq = freqi + freqihat;

        final double freqkappam1p1 = freq * (kappa - 1) + 1;

        return Math.exp(-t/N) * (pm * freqihat * Math.exp(mbetamutwotptau * freqkappam1p1) / (twobetamuN * freqkappam1p1 + 1) - freqj * ((1 - freq) * Math.exp(mbetamutwotptau) / (twobetamuN + 1) + freq)) / freq;

    }

    private double transversionH(final double t, final double N, final int j, final double tau, final double mu) {
        final double betamu = beta * mu;
        return frequencyModel.getFrequency(j) * Math.exp(-t/N) * (Math.exp(-betamu * (2*t + tau)) / (2 * betamu * N + 1) - 1);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object o, int i) {
        super.handleModelChangedEvent(model, o, i);
        betaKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int i, Variable.ChangeType changeType) {
        super.handleVariableChangedEvent(variable, i, changeType);
        betaKnown = false;
    }

    @Override
    protected void storeState() {
        // Nothing to do
    }

    @Override
    protected void restoreState() {
        super.restoreState();
        betaKnown = false;
    }

    @Override
    protected void acceptState() {
        // Nothing to do
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            return new HKYSkylineIntegrator((HKY) xo.getChild(HKY.class),
                    (PiecewisePopulationModel) xo.getChild(PiecewisePopulationModel.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        private final XMLSyntaxRule[] rules = {new ElementRule(HKY.class),
                new ElementRule(PiecewisePopulationModel.class)};

        @Override
        public String getParserDescription() {
            return "An integrator that supports the HKY and skyline models.";
        }

        @Override
        public Class getReturnType() {
            return HKYSkylineIntegrator.class;
        }

        @Override
        public String getParserName() {
            return "hkySkylineIntegrator";
        }
    };
}
