/*
 * TaxonSpecificSequenceErrorModel.java
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

import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.TaxonList;
import dr.evomodel.treelikelihood.SequenceErrorModel;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

import java.util.logging.Logger;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class TaxonSpecificSequenceErrorModel extends SequenceErrorModel {

    private final ErrorType errorType;
    private final Parameter baseErrorRateParameter;
    private final Parameter ageRelatedErrorRateParameter;
    private final Parameter indicatorParameter;

    public TaxonSpecificSequenceErrorModel(TaxonList includeTaxa, TaxonList excludeTaxa,
                                           ErrorType errorType, Parameter baseErrorRateParameter,
                                           Parameter ageRelatedErrorRateParameter,
                                           Parameter indicatorParameter) {
        super(includeTaxa, excludeTaxa, errorType, baseErrorRateParameter, ageRelatedErrorRateParameter, indicatorParameter);

        this.errorType = errorType;

        if (baseErrorRateParameter != null) {
            this.baseErrorRateParameter = baseErrorRateParameter;
            addVariable(this.baseErrorRateParameter);
        } else {
            this.baseErrorRateParameter = null;
        }

        if (ageRelatedErrorRateParameter != null) {
            this.ageRelatedErrorRateParameter = ageRelatedErrorRateParameter;
            addVariable(ageRelatedErrorRateParameter);
        } else {
            this.ageRelatedErrorRateParameter = null;
        }

        if (indicatorParameter != null) {
            this.indicatorParameter = indicatorParameter;
            addVariable(indicatorParameter);
        } else {
            this.indicatorParameter = null;
        }

    }

    @Override
    public void getTipPartials(int nodeIndex, double[] partials) {

        int[] states = this.states[nodeIndex];
        if (indicatorParameter == null || indicatorParameter.getParameterValue(nodeIndex) > 0.0) {

            double pUndamaged = 1.0;
            double pDamagedTS = 0.0;
            double pDamagedTV = 0.0;

            if (!excluded[nodeIndex]) {
                if (baseErrorRateParameter != null) {
                    pUndamaged = pUndamaged - baseErrorRateParameter.getParameterValue(nodeIndex);
                }

                if (ageRelatedErrorRateParameter != null) {
                    double rate = ageRelatedErrorRateParameter.getParameterValue(0);
                    double age = tree.getNodeHeight(tree.getExternalNode(nodeIndex));
                    pUndamaged *= Math.exp(-rate * age);
                }


                if (errorType == ErrorType.ALL_SUBSTITUTIONS) {
                    pDamagedTS = (1.0 - pUndamaged) / 3.0;
                    pDamagedTV = pDamagedTS;

                } else if (errorType == ErrorType.TRANSITIONS_ONLY) {
                    pDamagedTS = 1.0 - pUndamaged;
                    pDamagedTV = 0.0;
                } else {
                    throw new IllegalArgumentException("only TRANSITIONS_ONLY and ALL_SUBSTITUTIONS are supported");
                }

            }

            int k = 0;
            for (int j = 0; j < patternCount; j++) {
                switch (states[j]) {
                    case Nucleotides.A_STATE: // is an A
                        partials[k] = pUndamaged;
                        partials[k + 1] = pDamagedTV;
                        partials[k + 2] = pDamagedTS;
                        partials[k + 3] = pDamagedTV;
                        break;
                    case Nucleotides.C_STATE: // is an C
                        partials[k] = pDamagedTV;
                        partials[k + 1] = pUndamaged;
                        partials[k + 2] = pDamagedTV;
                        partials[k + 3] = pDamagedTS;
                        break;
                    case Nucleotides.G_STATE: // is an G
                        partials[k] = pDamagedTS;
                        partials[k + 1] = pDamagedTV;
                        partials[k + 2] = pUndamaged;
                        partials[k + 3] = pDamagedTV;
                        break;
                    case Nucleotides.UT_STATE: // is an T
                        partials[k] = pDamagedTV;
                        partials[k + 1] = pDamagedTS;
                        partials[k + 2] = pDamagedTV;
                        partials[k + 3] = pUndamaged;
                        break;
                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 1.0;
                }
                k += stateCount;
            }
        } else {
            int k = 0;
            for (int j = 0; j < patternCount; j++) {

                switch (states[j]) {
                    case Nucleotides.A_STATE: // is an A
                        partials[k] = 1.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                        break;
                    case Nucleotides.C_STATE: // is an C
                        partials[k] = 0.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                        break;
                    case Nucleotides.G_STATE: // is an G
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 0.0;
                        break;
                    case Nucleotides.UT_STATE: // is an T
                        partials[k] = 0.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 1.0;
                        break;
                    default: // is an ambiguity
                        partials[k] = 1.0;
                        partials[k + 1] = 1.0;
                        partials[k + 2] = 1.0;
                        partials[k + 3] = 1.0;
                }

                k += stateCount;
            }

        }
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String TAXON_SPECIFIC_SEQUENCE_ERROR_MODEL = "taxonSpecificSequenceErrorModel";
        public static final String BASE_ERROR_RATE = "baseErrorRate";
        public static final String AGE_RELATED_RATE = "ageRelatedErrorRate";
        public static final String INDICATORS = "indicators";

        public static final String EXCLUDE = "exclude";
        public static final String INCLUDE = "include";

        public static final String TYPE = "type";

        public String getParserName() {
            return TAXON_SPECIFIC_SEQUENCE_ERROR_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            SequenceErrorModel.ErrorType errorType = SequenceErrorModel.ErrorType.ALL_SUBSTITUTIONS;

            if (xo.hasAttribute(TYPE)) {
                if (xo.getStringAttribute(TYPE).equalsIgnoreCase("transitions")) {
                    errorType = SequenceErrorModel.ErrorType.TRANSITIONS_ONLY;
                } else if (!xo.getStringAttribute(TYPE).equalsIgnoreCase("all")) {
                    throw new XMLParseException("unrecognized option for attribute, 'type': " + xo.getStringAttribute(TYPE));
                }
            }

            Parameter baseDamageRateParameter = null;
            if (xo.hasChildNamed(BASE_ERROR_RATE)) {
                baseDamageRateParameter = (Parameter) xo.getElementFirstChild(BASE_ERROR_RATE);
            }

            Parameter ageRelatedRateParameter = null;
            if (xo.hasChildNamed(AGE_RELATED_RATE)) {
                ageRelatedRateParameter = (Parameter) xo.getElementFirstChild(AGE_RELATED_RATE);
            }

            if (baseDamageRateParameter == null && ageRelatedRateParameter == null) {
                throw new XMLParseException("You must specify one or other or both of " +
                        BASE_ERROR_RATE + " and " + AGE_RELATED_RATE + " parameters");
            }

            Parameter indicatorParameter = null;
            if (xo.hasChildNamed(INDICATORS)) {
                indicatorParameter = (Parameter)xo.getElementFirstChild(INDICATORS);
            }


            TaxonList includeTaxa = null;
            TaxonList excludeTaxa = null;

            if (xo.hasChildNamed(INCLUDE)) {
                includeTaxa = (TaxonList) xo.getElementFirstChild(INCLUDE);
            }

            if (xo.hasChildNamed(EXCLUDE)) {
                excludeTaxa = (TaxonList) xo.getElementFirstChild(EXCLUDE);
            }

            TaxonSpecificSequenceErrorModel aDNADamageModel = new TaxonSpecificSequenceErrorModel(includeTaxa, excludeTaxa,
                    errorType, baseDamageRateParameter, ageRelatedRateParameter, indicatorParameter);

            Logger.getLogger("dr.evomodel").info("Using sequence error model, assuming errors cause " +
                    (errorType == SequenceErrorModel.ErrorType.TRANSITIONS_ONLY ? "transitions only." : "any substitution."));

            return aDNADamageModel;
        }

        public String getParserDescription() {
            return "This element returns a model that allows for post-mortem DNA damage.";
        }

        public Class getReturnType() {
            return TaxonSpecificSequenceErrorModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(TYPE, true),
                new ElementRule(BASE_ERROR_RATE, Parameter.class, "The base error rate per site per sequence", true),
                new ElementRule(AGE_RELATED_RATE, Parameter.class, "The error rate per site per unit time", true),
                new ElementRule(INDICATORS, Parameter.class, "A binary indicator of whether the sequence has errors", true),
                new XORRule(
                        new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
                        new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to")
                        , true)
        };

    };

}
