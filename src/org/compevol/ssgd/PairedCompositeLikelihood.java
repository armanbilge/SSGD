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
import dr.evomodel.coalescent.PiecewisePopulationModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
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

    private final PatternList patterns;
    private final SiteModel siteModel;
    private final FrequencyModel frequencies;
    private final HKY hky;
    private final TipStatesModel tipStatesModel;
    private final PiecewisePopulationModel populationModel;
    private final int taxonCount;

    public PairedCompositeLikelihood(final PatternList patterns, final SiteModel siteModel, final TipStatesModel tipStatesModel, PiecewisePopulationModel populationModel) {
        super(new CompoundModel("PairedCompositeLikelihoodModel"));
        final CompoundModel model = (CompoundModel) getModel();
        this.patterns = patterns;
        this.siteModel = siteModel;
        model.addModel(siteModel);
        frequencies = siteModel.getFrequencyModel();
        hky = (HKY) siteModel.getSubstitutionModel();
        this.tipStatesModel = tipStatesModel;
        model.addModel(tipStatesModel);
        this.populationModel = populationModel;
        model.addModel(populationModel);
        taxonCount = patterns.getTaxonCount();

        // Hack involving a fake tree to set up the tip states model
        final ConstantPopulation constant = new ConstantPopulation(Units.Type.SUBSTITUTIONS);
        constant.setN0(1);
        final Tree fakeTree = new CoalescentSimulator().simulateTree(patterns, constant);
        tipStatesModel.setTree(fakeTree);
        for (int i = 0; i < fakeTree.getExternalNodeCount(); ++i) {
            final Taxon taxon = fakeTree.getNodeTaxon(fakeTree.getExternalNode(i));
            tipStatesModel.setStates(patterns, patterns.getTaxonIndex(taxon), i, taxon.getId());
        }
    }

    @Override
    protected double calculateLogLikelihood() {
        double logL = 0.0;

        for (int i = 0; i < patterns.getPatternCount(); ++i) {
            double patternLogL = 0.0;
            for (int j = 0; j < taxonCount; ++j)
                for (int k = j+1; k < taxonCount; ++k)
                    patternLogL += calculatePairLogLikelihood(j, k);
            logL += patternLogL * patterns.getPatternWeight(i);
        }

        return logL;
    }

    private double calculatePairLogLikelihood(final int i, final int j) {
        // TODO
        return 0.0;
    }

    private double beta() {
        final double kappa = hky.getKappa();
        final double A = frequencies.getFrequency(0);
        final double C = frequencies.getFrequency(1);
        final double G = frequencies.getFrequency(2);
        final double T = frequencies.getFrequency(3);
        return 1.0 / (2 * (A + G) * (C + T) + 2 * kappa * (A * G + C * T));
    }

    private double pAA() {
        // TODO
        return 0.0;
    }

    private double pAC() {
        // TODO
        return 0.0;
    }

    private double pAG() {
        // TODO
        return 0.0;
    }

    private double pAT() {
        // TODO
        return 0.0;
    }

    private double pCC() {
        // TODO
        return 0.0;
    }

    private double pCG() {
        // TODO
        return 0.0;
    }

    private double pCT() {
        // TODO
        return 0.0;
    }

    private double pGG() {
        // TODO
        return 0.0;
    }

    private double pGT() {
        // TODO
        return 0.0;
    }

    private double pTT() {
        // TODO
        return 0.0;
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new PairedCompositeLikelihood((PatternList) xo.getChild(PatternList.class),
                    (SiteModel) xo.getChild(SiteModel.class),
                    (TipStatesModel) xo.getChild(TipStatesModel.class),
                    (PiecewisePopulationModel) xo.getChild(PiecewisePopulationModel.class));
        }

        final XMLSyntaxRule[] rules = {new ElementRule(PatternList.class), new ElementRule(PatternList.class),
                new ElementRule(TipStatesModel.class), new ElementRule(PiecewisePopulationModel.class)};

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[0];
        }

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
