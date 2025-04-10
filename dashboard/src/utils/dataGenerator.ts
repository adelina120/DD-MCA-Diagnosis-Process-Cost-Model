import { Parameters, DataPoint } from '../types';

const utility = (pv: number, tat: number, alpha: number, lambda: number): number => {
  const beta = 1 - alpha;
  return alpha * (2 * pv - 1) + beta * Math.exp(-lambda * tat);
};

export const generateData = (params: Parameters, xAxis: keyof Parameters): DataPoint[] => {
  const {
    cmaCost,
    cmaYield,
    cmaPPV,
    cmaNPV,
    gpCost,
    gpYield,
    gpPPV,
    gpNPV,
    wesYield1Tier,
    wesYield2Tier,
    wesYield3Tier,
    wesCost,
    wesPPV,
    wesNPV,
    expertFee,
    alpha,
    lambda,
    aiPrecision,
    aiFDR,
    aiFOR,
    aiNPV
  } = params;

  const scenarios = [
    "Scenario 1 (CMA + GP)",
    "Scenario 2 (CMA + GP + WES)",
    "Scenario 3 (CMA + WES)",
    "Scenario 4 (WES alone)",
    "AI-delegation (r>r*)"
  ];

  // Generate values for the x-axis parameter
  const xAxisValues = Array.from({ length: 20 }, (_, i) => {
    const param = params[xAxis];
    if (typeof param === 'number') {
      // For parameters that are probabilities (yields, PPV, NPV, etc.)
      if (['cmaYield', 'gpYield', 'wesYield1Tier', 'wesYield2Tier', 'wesYield3Tier', 
           'cmaPPV', 'cmaNPV', 'gpPPV', 'gpNPV', 'wesPPV', 'wesNPV', 
           'aiPrecision', 'aiFDR', 'aiFOR', 'aiNPV'].includes(xAxis)) {
        return 0.1 + (i * 0.8) / 19; // Range from 0.1 to 0.9
      }
      // For cost parameters
      if (['cmaCost', 'gpCost', 'wesCost', 'expertFee'].includes(xAxis)) {
        return param * (0.5 + (i * 1.5) / 19); // Range from 50% to 200% of original value
      }
      // For alpha and lambda
      if (['alpha', 'lambda'].includes(xAxis)) {
        return 0.1 + (i * 0.8) / 19; // Range from 0.1 to 0.9
      }
    }
    return param;
  });

  return scenarios.flatMap(scenario => {
    return xAxisValues.map(xValue => {
      // Create a copy of parameters with the x-axis value updated
      const currentParams = { ...params, [xAxis]: xValue };

      // Scenario 1 (CMA + GP)
      const s1Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.gpCost + currentParams.expertFee);
      const s1Eff = currentParams.cmaYield * utility(currentParams.cmaPPV, 6, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.cmaYield) * (currentParams.gpYield * utility(currentParams.gpPPV, 14, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.gpYield) * utility(currentParams.gpNPV, 18, currentParams.alpha, currentParams.lambda));
      const s1EffCost = s1Cost / s1Eff;

      // Scenario 2 (CMA + GP + WES)
      const s2Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.gpCost + currentParams.expertFee + 
                    (1 - currentParams.gpYield) * (currentParams.wesCost + currentParams.expertFee));
      const s2Eff = currentParams.cmaYield * utility(currentParams.cmaPPV, 6, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.cmaYield) * currentParams.gpYield * utility(currentParams.gpPPV, 14, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.cmaYield) * (1 - currentParams.gpYield) * 
                    (currentParams.wesYield3Tier * utility(currentParams.wesPPV, 26, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.wesYield3Tier) * utility(currentParams.wesNPV, 26, currentParams.alpha, currentParams.lambda));
      const s2EffCost = s2Cost / s2Eff;

      // Scenario 3 (CMA + WES)
      const s3Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.wesCost + currentParams.expertFee);
      const s3Eff = currentParams.cmaYield * utility(currentParams.cmaPPV, 6, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.cmaYield) * (currentParams.wesYield2Tier * utility(currentParams.wesPPV, 18, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.wesYield2Tier) * utility(currentParams.wesNPV, 18, currentParams.alpha, currentParams.lambda));
      const s3EffCost = s3Cost / s3Eff;

      // Scenario 4 (WES alone)
      const s4Cost = currentParams.wesCost + currentParams.expertFee;
      const s4Eff = currentParams.wesYield1Tier * utility(currentParams.wesPPV, 12, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.wesYield1Tier) * utility(currentParams.wesNPV, 12, currentParams.alpha, currentParams.lambda);
      const s4EffCost = s4Cost / s4Eff;

      // AI-delegation mode
      const aiCost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.gpCost + (1 - currentParams.aiPrecision) * currentParams.wesCost);
      const aiEff = currentParams.aiPrecision * utility(currentParams.gpPPV, 10, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.aiPrecision) * (currentParams.wesYield3Tier * utility(currentParams.wesPPV, 18, currentParams.alpha, currentParams.lambda) + 
                    (1 - currentParams.wesYield3Tier) * utility(currentParams.wesNPV, 18, currentParams.alpha, currentParams.lambda));
      const aiEffCost = aiCost / aiEff;

      // Return the appropriate data point based on the scenario
      switch (scenario) {
        case "Scenario 1 (CMA + GP)":
          return {
            scenario,
            expectedCost: s1Cost,
            expectedUtility: s1Eff,
            effectiveCost: s1EffCost,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            alpha: currentParams.alpha,
            lambda: currentParams.lambda,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        case "Scenario 2 (CMA + GP + WES)":
          return {
            scenario,
            expectedCost: s2Cost,
            expectedUtility: s2Eff,
            effectiveCost: s2EffCost,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            alpha: currentParams.alpha,
            lambda: currentParams.lambda,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        case "Scenario 3 (CMA + WES)":
          return {
            scenario,
            expectedCost: s3Cost,
            expectedUtility: s3Eff,
            effectiveCost: s3EffCost,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            alpha: currentParams.alpha,
            lambda: currentParams.lambda,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        case "Scenario 4 (WES alone)":
          return {
            scenario,
            expectedCost: s4Cost,
            expectedUtility: s4Eff,
            effectiveCost: s4EffCost,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            alpha: currentParams.alpha,
            lambda: currentParams.lambda,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        case "AI-delegation (r>r*)":
          return {
            scenario,
            expectedCost: aiCost,
            expectedUtility: aiEff,
            effectiveCost: aiEffCost,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            alpha: currentParams.alpha,
            lambda: currentParams.lambda,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        default:
          throw new Error(`Unknown scenario: ${scenario}`);
      }
    });
  });
}; 