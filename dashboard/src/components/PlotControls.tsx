import React from 'react';
import { Box, FormControl, InputLabel, Select, MenuItem, Typography, SelectChangeEvent } from '@mui/material';
import { Parameters } from '../types';

interface PlotControlsProps {
  parameters: Parameters;
  onPlotConfigChange: (config: PlotConfig) => void;
}

export interface PlotConfig {
  scenarios: string[];
  yAxis: 'cost' | 'expectedUtility' | 'effectiveCost';
  xAxis: keyof Parameters;
  showAIComparison: boolean;
}

const availableScenarios = [
  "Scenario 1 (CMA + GP)",
  "Scenario 2 (CMA + GP + WES)",
  "Scenario 3 (CMA + WES)",
  "Scenario 4 (WES alone)",
  "AI-delegation (r>r*)"
];

const yAxisOptions = [
  { value: 'cost', label: 'Cost' },
  { value: 'expectedUtility', label: 'Expected Utility' },
  { value: 'effectiveCost', label: 'Effective Cost' }
];

const PlotControls: React.FC<PlotControlsProps> = ({ parameters, onPlotConfigChange }) => {
  const [config, setConfig] = React.useState<PlotConfig>({
    scenarios: availableScenarios,
    yAxis: 'effectiveCost',
    xAxis: 'alpha',
    showAIComparison: false
  });

  const handleScenarioChange = (event: SelectChangeEvent<string[]>) => {
    const newConfig = {
      ...config,
      scenarios: event.target.value as string[]
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  const handleYAxisChange = (event: SelectChangeEvent<'cost' | 'expectedUtility' | 'effectiveCost'>) => {
    const newConfig = {
      ...config,
      yAxis: event.target.value as 'cost' | 'expectedUtility' | 'effectiveCost'
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  const handleXAxisChange = (event: SelectChangeEvent<keyof Parameters>) => {
    const newConfig = {
      ...config,
      xAxis: event.target.value as keyof Parameters
    };
    setConfig(newConfig);
    onPlotConfigChange(newConfig);
  };

  // Get available x-axis parameters based on selected scenarios and y-axis
  const getAvailableXAxisParams = () => {
    const baseParams = [
      'cmaYield',
      'gpYield',
      'wesYield1Tier',
      'wesYield2Tier',
      'wesYield3Tier',
      'expertFee',
      'cmaCost',
      'gpCost',
      'wesCost'
    ];

    // Add AI-related parameters if AI scenario is selected
    if (config.scenarios.includes("AI-delegation (r>r*)")) {
      baseParams.push('aiPrecision', 'aiFDR', 'aiFOR', 'aiNPV');
    }

    // Add alpha and lambda only if y-axis is not cost
    if (config.yAxis !== 'cost') {
      baseParams.push('alpha', 'lambda');
    }

    return baseParams;
  };

  return (
    <Box sx={{ mb: 2 }}>
      <Typography gutterBottom>Plot Configuration</Typography>
      
      <FormControl fullWidth margin="normal">
        <InputLabel>Scenarios</InputLabel>
        <Select
          multiple
          value={config.scenarios}
          onChange={handleScenarioChange}
          label="Scenarios"
        >
          {availableScenarios.map(scenario => (
            <MenuItem key={scenario} value={scenario}>
              {scenario}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <FormControl fullWidth margin="normal">
        <InputLabel>Y-Axis</InputLabel>
        <Select
          value={config.yAxis}
          onChange={handleYAxisChange}
          label="Y-Axis"
        >
          {yAxisOptions.map(option => (
            <MenuItem key={option.value} value={option.value}>
              {option.label}
            </MenuItem>
          ))}
        </Select>
      </FormControl>

      <FormControl fullWidth margin="normal">
        <InputLabel>X-Axis</InputLabel>
        <Select
          value={config.xAxis}
          onChange={handleXAxisChange}
          label="X-Axis"
        >
          {getAvailableXAxisParams().map(param => (
            <MenuItem key={param} value={param}>
              {param}
            </MenuItem>
          ))}
        </Select>
      </FormControl>
    </Box>
  );
};

export default PlotControls; 