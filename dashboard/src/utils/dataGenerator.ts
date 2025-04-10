import { Parameters, DataPoint } from '../types';

const utility = (pv: number, tat: number, alpha: number, lambda: number): number => {
  const beta = 1 - alpha;
  return alpha * (2 * pv - 1) + beta * Math.exp(-lambda * tat);
};

export const generateData = (params: Parameters): DataPoint[] => {
  const {
    n,
    cmaYield,
    gpYield,
    wesYield1Tier,
    wesYield2Tier,
    wesYield3Tier,
    expertFee,
    cmaCost,
    gpCost,
    wesCost,
    cmaPPV,
    cmaNPV,
    gpPPV,
    gpNPV,
    wesPPV,
    wesNPV,
    aiPrecision,
    aiFDR,
    aiFOR,
    aiNPV,
    lambda
  } = params;

  // Generate multiple alpha values between 0.2 and 0.8
  const alphaValues = Array.from({ length: 20 }, (_, i) => 0.2 + (i * 0.6) / 19);

  const scenarios = [
    "Scenario 1 (CMA + GP)",
    "Scenario 2 (CMA + GP + WES)",
    "Scenario 3 (CMA + WES)",
    "Scenario 4 (WES alone)",
    "AI-delegation (r>r*)"
  ];

  return scenarios.flatMap(scenario => {
    return alphaValues.map(alpha => {
      // Scenario 1 (CMA + GP)
      const s1Cost = cmaCost + expertFee + (1 - cmaYield) * (gpCost + expertFee);
      const s1Eff = cmaYield * utility(cmaPPV, 6, alpha, lambda) + 
                    (1 - cmaYield) * (gpYield * utility(gpPPV, 14, alpha, lambda) + 
                    (1 - gpYield) * utility(gpNPV, 18, alpha, lambda));
      const s1EffCost = s1Cost / s1Eff;

      // Scenario 2 (CMA + GP + WES)
      const s2Cost = cmaCost + expertFee + (1 - cmaYield) * 
                    (gpCost + expertFee + (1 - gpYield) * (wesCost + expertFee));
      const s2Eff = cmaYield * utility(cmaPPV, 6, alpha, lambda) + 
                    (1 - cmaYield) * gpYield * utility(gpPPV, 14, alpha, lambda) + 
                    (1 - cmaYield) * (1 - gpYield) * 
                    (wesYield3Tier * utility(wesPPV, 26, alpha, lambda) + 
                    (1 - wesYield3Tier) * utility(wesNPV, 26, alpha, lambda));
      const s2EffCost = s2Cost / s2Eff;

      // Scenario 3 (CMA + WES)
      const s3Cost = cmaCost + expertFee + (1 - cmaYield) * (wesCost + expertFee);
      const s3Eff = cmaYield * utility(cmaPPV, 6, alpha, lambda) + 
                    (1 - cmaYield) * (wesYield2Tier * utility(wesPPV, 18, alpha, lambda) + 
                    (1 - wesYield2Tier) * utility(wesNPV, 18, alpha, lambda));
      const s3EffCost = s3Cost / s3Eff;

      // Scenario 4 (WES alone)
      const s4Cost = wesCost + expertFee;
      const s4Eff = wesYield1Tier * utility(wesPPV, 12, alpha, lambda) + 
                    (1 - wesYield1Tier) * utility(wesNPV, 12, alpha, lambda);
      const s4EffCost = s4Cost / s4Eff;

      // AI-delegation mode
      const aiCost = cmaCost + expertFee + (1 - cmaYield) * 
                    (gpCost + (1 - aiPrecision) * wesCost);
      const aiEff = aiPrecision * utility(gpPPV, 10, alpha, lambda) + 
                    (1 - aiPrecision) * (wesYield3Tier * utility(wesPPV, 18, alpha, lambda) + 
                    (1 - wesYield3Tier) * utility(wesNPV, 18, alpha, lambda));
      const aiEffCost = aiCost / aiEff;

      // Return the appropriate data point based on the scenario
      switch (scenario) {
        case "Scenario 1 (CMA + GP)":
          return {
            scenario,
            expectedCost: s1Cost,
            expectedUtility: s1Eff,
            effectiveCost: s1EffCost,
            cmaCost,
            gpCost,
            aiPerformance: aiPrecision,
            alphaValues: alpha
          };
        case "Scenario 2 (CMA + GP + WES)":
          return {
            scenario,
            expectedCost: s2Cost,
            expectedUtility: s2Eff,
            effectiveCost: s2EffCost,
            cmaCost,
            gpCost,
            aiPerformance: aiPrecision,
            alphaValues: alpha
          };
        case "Scenario 3 (CMA + WES)":
          return {
            scenario,
            expectedCost: s3Cost,
            expectedUtility: s3Eff,
            effectiveCost: s3EffCost,
            cmaCost,
            gpCost,
            aiPerformance: aiPrecision,
            alphaValues: alpha
          };
        case "Scenario 4 (WES alone)":
          return {
            scenario,
            expectedCost: s4Cost,
            expectedUtility: s4Eff,
            effectiveCost: s4EffCost,
            cmaCost,
            gpCost,
            aiPerformance: aiPrecision,
            alphaValues: alpha
          };
        case "AI-delegation (r>r*)":
          return {
            scenario,
            expectedCost: aiCost,
            expectedUtility: aiEff,
            effectiveCost: aiEffCost,
            cmaCost,
            gpCost,
            aiPerformance: aiPrecision,
            alphaValues: alpha
          };
        default:
          throw new Error(`Unknown scenario: ${scenario}`);
      }
    });
  });
}; 