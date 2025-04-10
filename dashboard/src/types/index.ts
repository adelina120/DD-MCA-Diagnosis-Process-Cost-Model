export interface Parameters {
  cmaCost: number;
  cmaYield: number;
  cmaPPV: number;
  cmaNPV: number;
  gpCost: number;
  gpYield: number;
  gpPPV: number;
  gpNPV: number;
  wesCost: number;
  wesYield: number;
  wesYield1Tier: number;
  wesYield2Tier: number;
  wesYield3Tier: number;
  wesPPV: number;
  wesNPV: number;
  expertFee: number;
  alpha: number;
  lambda: number;
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
  expectedUtility: number;
  effectiveCost: number;
  cmaCost: number;
  gpCost: number;
  aiPerformance: number;
  alphaValues: number;
  alpha: number;
  lambda: number;
  cmaYield: number;
  gpYield: number;
  wesYield1Tier: number;
  wesYield2Tier: number;
  wesYield3Tier: number;
  expertFee: number;
  wesCost: number;
  aiPrecision: number;
  aiFDR: number;
  aiFOR: number;
  aiNPV: number;
}

export interface PlotConfig {
  scenarios: string[];
  yAxis: 'cost' | 'expectedUtility' | 'effectiveCost';
  xAxis: keyof Parameters;
  showAIComparison: boolean;
}

export type PlotType = 'effectiveCost' | 'expectedUtility' | 'aiComparison'; 