import { Parameters, DataPoint } from '../types';

export const calculateResults = (parameters: Parameters): DataPoint[] => {
  const results: DataPoint[] = [];
  const alphaValues = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

  for (const alpha of alphaValues) {
    // Calculate for CMA
    const cmaCost = parameters.cmaCost;
    const cmaYield = parameters.cmaYield;
    const cmaUtility = cmaYield * (1 - alpha);

    // Calculate for GP
    const gpCost = parameters.gpCost;
    const gpYield = parameters.gpYield;
    const gpUtility = gpYield * (1 - alpha);

    // Calculate for WES
    const wesCost = parameters.wesCost;
    const wesYield = parameters.wesYield;
    const wesUtility = wesYield * (1 - alpha);

    // Calculate combined scenarios
    const cmaGpCost = cmaCost + gpCost;
    const cmaGpYield = 1 - (1 - cmaYield) * (1 - gpYield);
    const cmaGpUtility = cmaGpYield * (1 - alpha);

    const cmaWesCost = cmaCost + wesCost;
    const cmaWesYield = 1 - (1 - cmaYield) * (1 - wesYield);
    const cmaWesUtility = cmaWesYield * (1 - alpha);

    const cmaGpWesCost = cmaCost + gpCost + wesCost;
    const cmaGpWesYield = 1 - (1 - cmaYield) * (1 - gpYield) * (1 - wesYield);
    const cmaGpWesUtility = cmaGpWesYield * (1 - alpha);

    results.push({
      scenario: "Scenario 1 (CMA + GP)",
      expectedCost: cmaGpCost,
      expectedUtility: cmaGpUtility,
      effectiveCost: cmaGpCost / cmaGpUtility,
      cmaCost,
      gpCost,
      aiPerformance: parameters.aiPrecision,
      alphaValues: alpha,
      alpha: parameters.alpha,
      lambda: parameters.lambda,
      cmaYield: parameters.cmaYield,
      gpYield: parameters.gpYield,
      wesYield1Tier: parameters.wesYield1Tier,
      wesYield2Tier: parameters.wesYield2Tier,
      wesYield3Tier: parameters.wesYield3Tier,
      expertFee: parameters.expertFee,
      wesCost: parameters.wesCost,
      aiPrecision: parameters.aiPrecision,
      aiFDR: parameters.aiFDR,
      aiFOR: parameters.aiFOR,
      aiNPV: parameters.aiNPV
    });

    results.push({
      scenario: "Scenario 2 (CMA + GP + WES)",
      expectedCost: cmaGpWesCost,
      expectedUtility: cmaGpWesUtility,
      effectiveCost: cmaGpWesCost / cmaGpWesUtility,
      cmaCost,
      gpCost,
      aiPerformance: parameters.aiPrecision,
      alphaValues: alpha,
      alpha: parameters.alpha,
      lambda: parameters.lambda,
      cmaYield: parameters.cmaYield,
      gpYield: parameters.gpYield,
      wesYield1Tier: parameters.wesYield1Tier,
      wesYield2Tier: parameters.wesYield2Tier,
      wesYield3Tier: parameters.wesYield3Tier,
      expertFee: parameters.expertFee,
      wesCost: parameters.wesCost,
      aiPrecision: parameters.aiPrecision,
      aiFDR: parameters.aiFDR,
      aiFOR: parameters.aiFOR,
      aiNPV: parameters.aiNPV
    });

    results.push({
      scenario: "Scenario 3 (CMA + WES)",
      expectedCost: cmaWesCost,
      expectedUtility: cmaWesUtility,
      effectiveCost: cmaWesCost / cmaWesUtility,
      cmaCost,
      gpCost,
      aiPerformance: parameters.aiPrecision,
      alphaValues: alpha,
      alpha: parameters.alpha,
      lambda: parameters.lambda,
      cmaYield: parameters.cmaYield,
      gpYield: parameters.gpYield,
      wesYield1Tier: parameters.wesYield1Tier,
      wesYield2Tier: parameters.wesYield2Tier,
      wesYield3Tier: parameters.wesYield3Tier,
      expertFee: parameters.expertFee,
      wesCost: parameters.wesCost,
      aiPrecision: parameters.aiPrecision,
      aiFDR: parameters.aiFDR,
      aiFOR: parameters.aiFOR,
      aiNPV: parameters.aiNPV
    });

    results.push({
      scenario: "Scenario 4 (WES alone)",
      expectedCost: wesCost,
      expectedUtility: wesUtility,
      effectiveCost: wesCost / wesUtility,
      cmaCost,
      gpCost,
      aiPerformance: parameters.aiPrecision,
      alphaValues: alpha,
      alpha: parameters.alpha,
      lambda: parameters.lambda,
      cmaYield: parameters.cmaYield,
      gpYield: parameters.gpYield,
      wesYield1Tier: parameters.wesYield1Tier,
      wesYield2Tier: parameters.wesYield2Tier,
      wesYield3Tier: parameters.wesYield3Tier,
      expertFee: parameters.expertFee,
      wesCost: parameters.wesCost,
      aiPrecision: parameters.aiPrecision,
      aiFDR: parameters.aiFDR,
      aiFOR: parameters.aiFOR,
      aiNPV: parameters.aiNPV
    });

    results.push({
      scenario: "AI-delegation (r>r*)",
      expectedCost: cmaCost + gpCost * parameters.aiPrecision,
      expectedUtility: cmaUtility + gpUtility * parameters.aiPrecision,
      effectiveCost: (cmaCost + gpCost * parameters.aiPrecision) / (cmaUtility + gpUtility * parameters.aiPrecision),
      cmaCost,
      gpCost,
      aiPerformance: parameters.aiPrecision,
      alphaValues: alpha,
      alpha: parameters.alpha,
      lambda: parameters.lambda,
      cmaYield: parameters.cmaYield,
      gpYield: parameters.gpYield,
      wesYield1Tier: parameters.wesYield1Tier,
      wesYield2Tier: parameters.wesYield2Tier,
      wesYield3Tier: parameters.wesYield3Tier,
      expertFee: parameters.expertFee,
      wesCost: parameters.wesCost,
      aiPrecision: parameters.aiPrecision,
      aiFDR: parameters.aiFDR,
      aiFOR: parameters.aiFOR,
      aiNPV: parameters.aiNPV
    });
  }

  return results;
}; 