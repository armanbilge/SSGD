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

import dr.evolution.alignment.Patterns;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
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

import java.util.Arrays;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class PairedCompositeLikelihood extends Likelihood.Abstract {

    private PairedPatterns patterns;
    private final SiteModel siteModel;
    private final Integrator integrator;
    private final TipStatesModel tipStatesModel;
    private final TaxonList taxa;

    public PairedCompositeLikelihood(final PairedPatterns patterns, final SiteModel siteModel, final Integrator integrator, final TipStatesModel tipStatesModel) {
        super(new CompoundModel("PairedCompositeLikelihoodModel"));
        final CompoundModel model = (CompoundModel) getModel();
        this.patterns = patterns;
        this.siteModel = siteModel;
        model.addModel(siteModel);
        this.integrator = integrator;
        model.addModel(integrator);
        this.tipStatesModel = tipStatesModel;
        model.addModel(tipStatesModel);
        taxa = patterns.getTaxa();

        // Hack involving fake patterns and a fake tree to set up the tip states model
        final Patterns fakePatterns = new Patterns(patterns.getDataType(), taxa);
        for (int i = 0; i < patterns.getDataType().getStateCount(); ++i) {
            final int[] pattern = new int[taxa.getTaxonCount()];
            Arrays.fill(pattern, i);
            fakePatterns.addPattern(pattern);
        }
        SimpleNode root = new SimpleNode();
        root.setTaxon(taxa.getTaxon(0));
        root.setHeight(root.getTaxon().getHeight());
        for (int i = 1; i < taxa.getTaxonCount(); ++i) {
            final SimpleNode child = new SimpleNode();
            child.setTaxon(taxa.getTaxon(i));
            child.setHeight(child.getTaxon().getHeight());
            final SimpleNode newRoot = new SimpleNode();
            newRoot.addChild(root);
            newRoot.addChild(child);
            root = newRoot;
        }
        final Tree fakeTree = new SimpleTree(root);
        tipStatesModel.setTree(fakeTree);
        for (int i = 0; i < taxa.getTaxonCount(); ++i)
            tipStatesModel.setStates(fakePatterns, i, i, taxa.getTaxon(i).getId());
    }

    @Override
    protected double calculateLogLikelihood() {

        final int taxonCount = taxa.getTaxonCount();
        final int stateCount = patterns.getDataType().getStateCount();

        final double[][] partials = new double[taxonCount][stateCount * stateCount];
        for (int i = 0; i < taxonCount; ++i) {
            tipStatesModel.getTipPartials(i, partials[i]);
        }

        double logL = 0.0;

        for (int x = 0; x < taxonCount; ++x) {

            for (int y = x+1; y < taxonCount; ++y) {

                for (int i = 0; i < stateCount; ++i) {

                    for (int j = 0; j < stateCount; ++j) {

                        final Taxon a = taxa.getTaxon(x);
                        final double[] aPartial = new double[stateCount];
                        System.arraycopy(partials[i], stateCount * i, aPartial, 0, stateCount);

                        final Taxon b = taxa.getTaxon(y);
                        final double[] bPartial = new double[stateCount];
                        System.arraycopy(partials[i], stateCount * j, bPartial, 0, stateCount);

                        logL += patterns.getPatternWeight(a, i, b, j) * pairLogLikelihood(a, aPartial, b, bPartial);

                    }

                }

            }

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
                    categoryL += frequencies.getFrequency(i) * aPartial[i] * bPartial[j] * integrator.integratedProbability(i, a.getHeight(), j, b.getHeight(), mu);
                }
            }

            L += categoryL * siteModel.getProportionForCategory(c);
        }

        return Math.log(L);

    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new PairedCompositeLikelihood((PairedPatterns) xo.getChild(PairedPatterns.class),
                    (SiteModel) xo.getChild(SiteModel.class),
                    (Integrator) xo.getChild(Integrator.class),
                    (TipStatesModel) xo.getChild(TipStatesModel.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(PairedPatterns.class), new ElementRule(SiteModel.class),
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
