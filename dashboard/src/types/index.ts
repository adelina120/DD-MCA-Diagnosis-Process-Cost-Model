export interface Parameters {
  cmaCost: number;
  cmaYield: number;
  cmaPPV: number;
  cmaNPV: number;
  cmaTAT: number;
  gpCost: number;
  gpYield: number;
  gpPPV: number;
  gpNPV: number;
  gpTAT: number;
  wesCost: number;
  wesYield: number;
  wesYield1Tier: number;
  wesYield2Tier: number;
  wesYield3Tier: number;
  wesPPV: number;
  wesNPV: number;
  wesTAT: number;
  expertFee: number;
  expertTAT: number;
  uTP: number;
  uFP: number;
  uTN: number;
  uFN: number;
  uInitial: number;
  numberOfYears: number;
  aiPrecision: number;
  aiFDR: number;
  aiFOR: number;
  aiNPV: number;
}

export interface ScenarioResult {
  cost: number;
  yield: number;
  utility: number;
}

export interface DataPoint {
  scenario: string;
  expectedCost: number;
  expectedQALY: number;
  costPerQALY: number;
  cmaCost: number;
  gpCost: number;
  aiPerformance: number;

  cmaYield: number;
  gpYield: number;
  wesYield1Tier: number;
  wesYield2Tier: number;
  wesYield3Tier: number;
  cmaPPV: number;
  cmaNPV: number;
  gpPPV: number;
  gpNPV: number;
  wesPPV: number;
  wesNPV: number;
  expertFee: number;
  wesCost: number;
  uTP: number;
  uFP: number;
  uTN: number;
  uFN: number;
  uInitial: number;
  numberOfYears: number;
  aiPrecision: number;
  aiFDR: number;
  aiFOR: number;
  aiNPV: number;
}

export interface PlotConfig {
  scenarios: string[];
  yAxis: 'cost' | 'expectedQALY' | 'costPerQALY';
  xAxis: keyof Parameters;
  showAIComparison: boolean;
}

export type PlotType = 'costPerQALY' | 'expectedQALY' | 'aiComparison'; 