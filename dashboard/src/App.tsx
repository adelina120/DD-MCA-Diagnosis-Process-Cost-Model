import React, { useState, useEffect } from 'react';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import { Box, Container, Grid, Paper, Typography } from '@mui/material';
import ParameterControls from './components/ParameterControls';
import Plots from './components/Plots';
import PlotControls from './components/PlotControls';
import { Parameters, PlotConfig, DataPoint } from './types';
import { generateData } from './utils/dataGenerator';

const theme = createTheme({
  palette: {
    mode: 'light',
  },
});

const initialParameters: Parameters = {
  cmaCost: 1000,
  cmaYield: 0.8,
  cmaPPV: 0.9,
  cmaNPV: 0.95,
  cmaTAT: 2,
  gpCost: 2000,
  gpYield: 0.7,
  gpPPV: 0.85,
  gpNPV: 0.9,
  gpTAT: 4, 
  wesCost: 3000,
  wesYield: 0.9,
  wesYield1Tier: 0.6,
  wesYield2Tier: 0.7,
  wesYield3Tier: 0.8,
  wesPPV: 0.95,
  wesNPV: 0.98,
  wesTAT: 8.0, // in weeks
  expertFee: 500,
  expertTAT: 4, // in weeks
  uTP: 0.6,
  uFP: 0.5,
  uTN: 1.0,
  uFN: -0.1,
  uInitial: 0,
  numberOfYears: 1.0,
  aiPrecision: 0.85,
  aiFDR: 0.1,
  aiFOR: 0.05,
  aiNPV: 0.95,
  cmaTAT: 2,        // Add these missing properties
  gpTAT: 4,
  wesTAT: 8,
  expertTAT: 4,
};

const initialPlotConfig: PlotConfig = {
  scenarios: [
    "Expert-alone: Scenario 1 (CMA + GP)",
    "Expert-alone: Scenario 2 (CMA + GP + WES)",
    "Expert-alone: Scenario 3 (CMA + WES)",
    "Expert-alone: Scenario 4 (WES alone)",
    "AI-delegation: r>r* (CMA + GP + WES)"
  ],
  yAxis: 'costPerQALY',
  xAxis: 'numberOfYears',
  showAIComparison: false
};

const App: React.FC = () => {
  const [parameters, setParameters] = useState<Parameters>(initialParameters);
  const [plotConfig, setPlotConfig] = useState<PlotConfig>(initialPlotConfig);
  const [results, setResults] = useState<DataPoint[]>([]);

  useEffect(() => {
    const newResults = generateData(parameters, plotConfig.xAxis);
    setResults(newResults);
  }, [parameters, plotConfig.xAxis]);

  const handleParametersChange = (newParameters: Parameters) => {
    setParameters(newParameters);
  };

  const handlePlotConfigChange = (newPlotConfig: PlotConfig) => {
    setPlotConfig(newPlotConfig);
  };

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <Container maxWidth="xl">
        <Box sx={{ my: 4 }}>
          <Typography variant="h4" component="h1" gutterBottom>
            MCA Diagnosis Process Cost Model Dashboard
          </Typography>
          <Grid container spacing={3}>
            <Grid item xs={12} md={4}>
              <Paper sx={{ p: 2 }}>
                <ParameterControls
                  parameters={parameters}
                  onParametersChange={handleParametersChange}
                  xAxis={plotConfig.xAxis}
                />
              </Paper>
            </Grid>
            <Grid item xs={12} md={8}>
              <Paper sx={{ p: 2 }}>
                <PlotControls
                  parameters={parameters}
                  onPlotConfigChange={handlePlotConfigChange}
                />
                <Box sx={{ mt: 2 }}>
                  <Plots
                    data={results}
                    plotConfig={plotConfig}
                  />
                </Box>
              </Paper>
            </Grid>
          </Grid>
        </Box>
      </Container>
    </ThemeProvider>
  );
};

export default App; 