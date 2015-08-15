/*
 * PairedCompositeLikelihood.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.coalescent.CoalescentSimulator;
import dr.evolution.coalescent.ConstantPopulation;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.treelikelihood.TipStatesModel;
import dr.inference.model.CompoundModel;
import dr.inference.model.Likelihood;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class PairedCompositeLikelihood extends Likelihood.Abstract {

    private PatternList patterns;
    private final SiteModel siteModel;
    private final Integrator integrator;
    private final TipStatesModel tipStatesModel;

    public PairedCompositeLikelihood(final PatternList patterns, final SiteModel siteModel, final Integrator integrator, final TipStatesModel tipStatesModel) {
        super(new CompoundModel("PairedCompositeLikelihoodModel"));
        final CompoundModel model = (CompoundModel) getModel();
        this.patterns = patterns;
        this.siteModel = siteModel;
        model.addModel(siteModel);
        this.integrator = integrator;
        model.addModel(integrator);
        this.tipStatesModel = tipStatesModel;
        model.addModel(tipStatesModel);

        // Hack involving a fake tree to set up the tip states model
        final ConstantPopulation constant = new ConstantPopulation(Units.Type.SUBSTITUTIONS);
        constant.setN0(1);
        final Tree fakeTree = new CoalescentSimulator().simulateTree(patterns, constant);
        tipStatesModel.setTree(fakeTree);
        for (int i = 0; i < patterns.getTaxonCount(); ++i)
            tipStatesModel.setStates(patterns, i, i, patterns.getTaxon(i).getId());
    }

    public void setPatterns(final PatternList patterns) {
        this.patterns = patterns;
        for (int i = 0; i < patterns.getTaxonCount(); ++i)
            tipStatesModel.setStates(patterns, i, i, patterns.getTaxon(i).getId());
    }

    @Override
    protected double calculateLogLikelihood() {

        double logL = 0.0;

        final int taxonCount = patterns.getTaxonCount();
        final int patternCount = patterns.getPatternCount();
        final int stateCount = siteModel.getSubstitutionModel().getDataType().getStateCount();

        final double[][] partials = new double[taxonCount][patternCount * stateCount];
        for (int i = 0; i < taxonCount; ++i) {
            tipStatesModel.getTipPartials(i, partials[i]);
        }

        int l = 0;
        for (int i = 0; i < patternCount; ++i) {

            double patternLogL = 0.0;

            for (int j = 0; j < taxonCount; ++j) {

                for (int k = j+1; k < taxonCount; ++k) {

                    final Taxon a = patterns.getTaxon(j);
                    final double[] aPartial = new double[stateCount];
                    System.arraycopy(partials[j], l, aPartial, 0, stateCount);
                    final Taxon b = patterns.getTaxon(k);
                    final double[] bPartial = new double[stateCount];
                    System.arraycopy(partials[k], l, bPartial, 0, stateCount);

                    patternLogL += pairLogLikelihood(a, aPartial, b, bPartial);

                }

            }

            logL += patternLogL * patterns.getPatternWeight(i);
            l += stateCount;
        }

        return logL;
    }

    private double pairLogLikelihood(final Taxon a, final double[] aPartial, final Taxon b, final double[] bPartial) {

        double L = 0.0;

        for (int c = 0; c < siteModel.getCategoryCount(); ++c) {

            final double mu = siteModel.getRateForCategory(c);

            double categoryL = 0.0;

            final int stateCount = siteModel.getSubstitutionModel().getDataType().getStateCount();
            final FrequencyModel frequencies = siteModel.getFrequencyModel();

            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    categoryL += frequencies.getFrequency(i)* aPartial[i] * bPartial[j] * integrator.integratedProbability(i, a.getHeight(), j, b.getHeight(), mu);
                }
            }

            L += categoryL * siteModel.getProportionForCategory(c);
        }

        return Math.log(L);

    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new PairedCompositeLikelihood((PatternList) xo.getChild(PatternList.class),
                    (SiteModel) xo.getChild(SiteModel.class),
                    (Integrator) xo.getChild(Integrator.class),
                    (TipStatesModel) xo.getChild(TipStatesModel.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(PatternList.class), new ElementRule(SiteModel.class),
                new ElementRule(Integrator.class), new ElementRule(TipStatesModel.class)};


        @Override
        public String getParserDescription() {
            return "Likelihood of a sequence evolution model using a sequence-pair approximation.";
        }

        @Override
        public Class getReturnType() {
            return PairedCompositeLikelihood.class;
        }

        @Override
        public String getParserName() {
            return "pairedCompositeLikelihood";
        }
    };

}
