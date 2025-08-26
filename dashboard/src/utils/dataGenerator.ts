import { Parameters, DataPoint } from '../types';

const utility = (pv: number, uT: number, uF: number): number => {
  return (pv*uT + (1-pv)*uF);
};

const qaly = (x: number, tat: number, uS0: number, eu: number): number=> {
  return(tat*uS0 + (x-tat)*eu);
}

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
    uTP,
    uFP,
    uTN,
    uFN,
    uInitial,
    numberOfYears,
    aiPrecision,
    aiFDR,
    aiFOR,
    aiNPV
  } = params;

  const scenarios = [
    "Expert-alone: Scenario 1 (CMA + GP)",
    "Expert-alone: Scenario 2 (CMA + GP + WES)",
    "Expert-alone: Scenario 3 (CMA + WES)",
    "Expert-alone: Scenario 4 (WES alone)",
    "AI-delegation: r>r* (CMA + GP + WES)"
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
      // For utility parameters
      if (['uTP', 'uFP', 'uTN', 'uFN', 'uInitial'].includes(xAxis)) {
        return -1 + (i * 2) / 19; // Range from -1 to 1
      }
    }
    return param;
  });

  return scenarios.flatMap(scenario => {
    return xAxisValues.map(xValue => {
      // Create a copy of parameters with the x-axis value updated
      const currentParams = { ...params, [xAxis]: xValue };

      // Calculate expected utilites for each ending state
      const uCMApos = utility(currentParams.cmaPPV, currentParams.uTP, currentParams.uFP);
      const uGPpos = utility(currentParams.gpPPV, currentParams.uTP, currentParams.uFP);
      const uWESpos = utility(currentParams.wesPPV, currentParams.uTP, currentParams.uFP);
      const uGPneg = utility(currentParams.gpNPV, currentParams.uTN,  currentParams.uFN);
      const uWESneg = utility(currentParams.wesNPV, currentParams.uTN, currentParams.uFN);

      // Turnaround times for each possible path
      const cmaTAT = (6*7)/365; // 6 weeks
      const cmaGpTAT = (14*7)/365; // 14 weeks
      const cmaWesTAT = (18*7)/365; // 18 weeks
      const cmaGpWesTAT = (26*7)/365; // 18 weeks
      const wesTAT = (12*7)/365; // 12 weeks
      const aiTAT1 = (10*7)/365; // 10 weeks
      const aiTAT2 = (18*7)/365; // 18 weeks

      // Scenario 1 (CMA + GP)
      const s1Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield)*(currentParams.gpCost + currentParams.expertFee);
      const s1QALY = currentParams.cmaYield*qaly(numberOfYears,cmaTAT,uInitial,uCMApos) + 
                    (1 - currentParams.cmaYield) * (currentParams.gpYield * qaly(numberOfYears,cmaGpTAT,uInitial,uGPpos) +
                    (1 - currentParams.gpYield) * qaly(numberOfYears,cmaGpTAT,uInitial,uGPneg));
      const s1Cpq = s1Cost / s1QALY;
      
      // Scenario 2 (CMA + GP + WES)
      const s2Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.gpCost + currentParams.expertFee + 
                    (1 - currentParams.gpYield) * (currentParams.wesCost + currentParams.expertFee));
      const s2QALY = cmaYield * qaly(numberOfYears, cmaTAT, uInitial, uCMApos)
              + (1-cmaYield) * gpYield * qaly(numberOfYears, cmaGpTAT, uInitial, uGPpos)
              + (1-cmaYield) * (1-gpYield) * wesYield3Tier * qaly(numberOfYears, cmaGpWesTAT, uInitial, uWESpos)
              + (1-cmaYield) * (1-gpYield) * (1-wesYield3Tier) * qaly(numberOfYears, cmaGpWesTAT, uInitial, uWESneg);
      const s2Cpq = s2Cost / s2QALY;

      // Scenario 3 (CMA + WES)
      const s3Cost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.wesCost + currentParams.expertFee);
      const s3QALY = currentParams.cmaYield * qaly(numberOfYears, cmaTAT, uInitial, uCMApos) + 
                    (1 - currentParams.cmaYield) * (currentParams.wesYield2Tier * qaly(numberOfYears, cmaWesTAT, uInitial, uWESpos) + 
                    (1 - currentParams.wesYield2Tier) * qaly(numberOfYears, cmaWesTAT, uInitial, uWESneg));
      const s3Cpq = s3Cost / s3QALY;

      // Scenario 4 (WES alone)
      const s4Cost = currentParams.wesCost + currentParams.expertFee;
      const s4QALY = currentParams.wesYield1Tier * qaly(numberOfYears, wesTAT, uInitial, uWESpos) +
                    (1 - currentParams.wesYield1Tier) * qaly(numberOfYears, wesTAT, uInitial, uWESneg);
      const s4Cpq = s4Cost / s4QALY;

      // AI-delegation mode
      const aiCost = currentParams.cmaCost + currentParams.expertFee + 
                    (1 - currentParams.cmaYield) * (currentParams.gpCost + (1 - currentParams.aiPrecision) * currentParams.wesCost);
      const aiQALY = currentParams.cmaYield * qaly(numberOfYears, cmaTAT, uInitial, uCMApos) + 
                    (1 - currentParams.cmaYield) * (currentParams.aiPrecision * qaly(numberOfYears, aiTAT1, uInitial, uGPpos) + 
                    (1 - currentParams.aiPrecision) * (currentParams.wesYield3Tier * qaly(numberOfYears, aiTAT2, uInitial, uWESpos) + 
                    (1 - currentParams.wesYield3Tier) * qaly(numberOfYears, aiTAT2, uInitial, uWESneg)));
      const aiCpq = aiCost / aiQALY;

      // Return the appropriate data point based on the scenario
      switch (scenario) {
        case "Scenario 1 (CMA + GP)":
          return {
            scenario,
            expectedCost: s1Cost,
            expectedQALY: s1QALY,
            costPerQALY: s1Cpq,
            cmaCost: currentParams.cmaCost,
            gpCost: currentParams.gpCost,
            aiPerformance: currentParams.aiPrecision,
            alphaValues: xValue,
            cmaYield: currentParams.cmaYield,
            gpYield: currentParams.gpYield,
            wesYield1Tier: currentParams.wesYield1Tier,
            wesYield2Tier: currentParams.wesYield2Tier,
            wesYield3Tier: currentParams.wesYield3Tier,
            cmaPPV: currentParams.cmaPPV,
            cmaNPV: currentParams.cmaNPV,
            gpPPV: currentParams.gpPPV,
            gpNPV: currentParams.gpNPV,
            wesPPV: currentParams.wesPPV,
            wesNPV: currentParams.wesNPV,
            expertFee: currentParams.expertFee,
            wesCost: currentParams.wesCost,
            uTP: currentParams.uTP,
            uFP: currentParams.uFP,
            uTN: currentParams.uTN,
            uFN: currentParams.uFN,
            uInitial: currentParams.uInitial,
            numberOfYears: currentParams.numberOfYears,
            aiPrecision: currentParams.aiPrecision,
            aiFDR: currentParams.aiFDR,
            aiFOR: currentParams.aiFOR,
            aiNPV: currentParams.aiNPV
          };
        case "Scenario 2 (CMA + GP + WES)":
          return {
            scenario,
            expectedCost: s2Cost,
            expectedEffectiveness: s2Eff,
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
            cmaPPV: currentParams.cmaPPV,
            cmaNPV: currentParams.cmaNPV,
            gpPPV: currentParams.gpPPV,
            gpNPV: currentParams.gpNPV,
            wesPPV: currentParams.wesPPV,
            wesNPV: currentParams.wesNPV,
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
            expectedEffectiveness: s3Eff,
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
            cmaPPV: currentParams.cmaPPV,
            cmaNPV: currentParams.cmaNPV,
            gpPPV: currentParams.gpPPV,
            gpNPV: currentParams.gpNPV,
            wesPPV: currentParams.wesPPV,
            wesNPV: currentParams.wesNPV,
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
            expectedEffectiveness: s4Eff,
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
            cmaPPV: currentParams.cmaPPV,
            cmaNPV: currentParams.cmaNPV,
            gpPPV: currentParams.gpPPV,
            gpNPV: currentParams.gpNPV,
            wesPPV: currentParams.wesPPV,
            wesNPV: currentParams.wesNPV,
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
            expectedEffectiveness: aiEff,
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
            cmaPPV: currentParams.cmaPPV,
            cmaNPV: currentParams.cmaNPV,
            gpPPV: currentParams.gpPPV,
            gpNPV: currentParams.gpNPV,
            wesPPV: currentParams.wesPPV,
            wesNPV: currentParams.wesNPV,
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
