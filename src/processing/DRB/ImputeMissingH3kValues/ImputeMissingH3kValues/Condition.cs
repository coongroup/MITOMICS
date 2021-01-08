using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Distributions;

namespace ImputeMissingH3kValues
{
    public class Condition
    {
        public static List<string> biomolecules;
        public string conditionName;
        public int imputationCutoff;
        public List<Replicate> replicates;
        public List<LfqValues> valuesToImputeFrom;
        public List<LfqValues> valuesToImpute;

        public Condition(string conditionName)
        {
            this.conditionName = conditionName;
            replicates = new List<Replicate>();
        }

        public void RunImputationAlgorithm()
        {
            var bestCv = new TargetCv();

            if (biomolecules == null)
            {
                biomolecules = GetDistinctBiomolecules();
            }

            var biomoleculeDictionary = CreateLfqDictionary(biomolecules);
            var nonImputedValues = new List<LfqValues>();
            var imputedValues = new List<LfqValues>();


            ExtractDataForImputation(out nonImputedValues, out imputedValues);

            // Data viewer algorithm
            // for (var i = 1; i <= 1000; i++)
            // use this one to speed up algorithm for a specific imputation cutoff
            for (var i = 100; i <= 101; i++)
            {
                var calculatedCv = Impute(i, biomoleculeDictionary, nonImputedValues, imputedValues);

                if (calculatedCv < bestCv.cv)
                {
                    bestCv.cv = calculatedCv;
                    bestCv.cutoff = i;
                }

                if (i == 1 || i == 90 || i == 1000)
                {
                    var t = "";
                }
            }
            

            // reimpute for best cv
            Impute(bestCv.cutoff, biomoleculeDictionary, nonImputedValues, imputedValues);
            // impute for a cutoff of 10%
            Impute(100, biomoleculeDictionary, nonImputedValues, imputedValues);
            this.imputationCutoff = bestCv.cutoff;
        }

        public double Impute(int cutoff, Dictionary<string, List<LfqValues>> biomoleculeDictionary, List<LfqValues> nonImputedValues, List<LfqValues> imputedValues)
        {
            var distribution = ExtractDistributionValues(nonImputedValues, cutoff);

            foreach (var value in imputedValues)
            {
                var sampled = distribution.Sample();

                while (sampled < 0)
                {
                    sampled = distribution.Sample();
                }

                value.lfq = sampled;
            }

            return CalculateAverageCV(biomoleculeDictionary);
        }

        private void ExtractDataForImputation(out List<LfqValues> nonImputedValues, out List<LfqValues> imputedValues)
        {
            // grab non-imputed values
            nonImputedValues = new List<LfqValues>();
            foreach (var replicate in replicates)
            {
                nonImputedValues = nonImputedValues.Concat(replicate.lfqValues.Where(lfqValue => !lfqValue.isImputed).ToList()).ToList();
            }
            nonImputedValues.Sort((a, b) => a.lfq.CompareTo(b.lfq));

            imputedValues = new List<LfqValues>();
            foreach (var replicate in replicates)
            {
                imputedValues = imputedValues.Concat(replicate.lfqValues.Where(lfqValue => lfqValue.isImputed).ToList()).ToList();
            }
        }

        private Normal ExtractDistributionValues(List<LfqValues> nonImputedValues, int cutoff)
        {
            var distributionValues = nonImputedValues.Select(lfq => lfq.lfq).ToList();

            var numToTake = (int)Math.Floor(distributionValues.Count * (cutoff / 1000.0));

            var values = distributionValues.GetRange(0, numToTake);

            var mean = values.Average();
            var stDev = CalculateStdDev(values);

            var distribution = new Normal(mean, stDev, Global.Random);

            return distribution;
        }

        private double CalculateMean(List<double> values)
        {
            var tempValues = new List<double>();

            foreach (var value in values)
            {
                tempValues.Add(Math.Pow(2, value));
            }

            double ret = 0;

            if (tempValues.Count() > 0)
            {
                //Compute the Average      
                ret = tempValues.Average();
            }
            return ret;
        }

        private double CalculateStdDev(List<double> values)
        {
            double ret = 0;

            if (values.Count() > 0)
            {
                //Compute the Average      
                double avg = values.Average();
                //Perform the Sum of (value-avg)_2_2      
                double sum = values.Sum(d => Math.Pow(d - avg, 2));
                //Put it all together      
                ret = Math.Sqrt((sum) / (values.Count() - 1));
            }
            return ret;
        }

        private double CalculateStdDevLog2Transformed(List<double> values)
        {
            var tempValues = new List<double>();

            foreach (var value in values)
            {
                tempValues.Add(Math.Pow(2, value));
            }

            double ret = 0;

            if (tempValues.Count() > 0)
            {
                //Compute the Average      
                double avg = tempValues.Average();
                //Perform the Sum of (value-avg)_2_2      
                double sum = tempValues.Sum(d => Math.Pow(d - avg, 2));
                //Put it all together      
                ret = Math.Sqrt((sum) / (tempValues.Count() - 1));
            }
            return ret;
        }

        private double CalculateAverageCV(Dictionary<string, List<LfqValues>> biomoleculeDictionary)
        {
            var cvList = new List<double>();
            var biomolecules = biomoleculeDictionary.Keys;

            foreach (var biomolecule in biomolecules)
            {
                var biomoleculeLfqValues = biomoleculeDictionary[biomolecule];

                cvList.Add(CalculateCV(biomoleculeLfqValues));
            }

            return cvList.Average();
        }

        private double CalculateMeanCV(Dictionary<string, List<LfqValues>> biomoleculeDictionary)
        {
            var cvList = new List<double>();
            var biomolecules = biomoleculeDictionary.Keys;

            foreach (var biomolecule in biomolecules)
            {
                var biomoleculeLfqValues = biomoleculeDictionary[biomolecule];

                cvList.Add(CalculateCV(biomoleculeLfqValues));
            }

            return 0;
        }

        private double CalculateCV(List<LfqValues> biomoleculeLfqValues)
        {
            var lfqValues = biomoleculeLfqValues.Select(x => x.lfq).ToList();

            var mean = 0.0;

            foreach (var x in lfqValues)
            {
                mean += Math.Pow(2, x);
            }
            mean = mean / lfqValues.Count;

            var stDev = CalculateStdDevLog2Transformed(lfqValues);

            return stDev / mean * 100;
        }

        private List<string> GetDistinctBiomolecules()
        {
            var returnList = new List<string>();

            var allValues = new List<LfqValues>();

            foreach (var replicate in replicates)
            {
                allValues = allValues.Concat(replicate.lfqValues).ToList();
            }

            returnList = allValues.Select(lfqvalue => lfqvalue.bioMolecule).Distinct().ToList();

            return returnList;
        }

        private Dictionary<string, List<LfqValues>> CreateLfqDictionary(List<string> biomolecules)
        {
            var returnDictionary = new Dictionary<string, List<LfqValues>>();
            var allValues = new List<LfqValues>();

            foreach (var replicate in replicates)
            {
                allValues = allValues.Concat(replicate.lfqValues).ToList();
            }

            foreach (var biomolecule in biomolecules)
            {
                var lfqValues = allValues.Where(lfqValue => lfqValue.bioMolecule.Equals(biomolecule)).ToList();
                returnDictionary.Add(biomolecule, lfqValues);
            }

            return returnDictionary;
        }
    }
}
