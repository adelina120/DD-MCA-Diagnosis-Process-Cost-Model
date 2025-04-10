import React, { useState } from 'react';
import { ThemeProvider, createTheme } from '@mui/material/styles';
import CssBaseline from '@mui/material/CssBaseline';
import { Box, Container, Grid, Paper, Typography } from '@mui/material';
import ParameterControls from './components/ParameterControls';
import Plots from './components/Plots';
import { Parameters } from './types';
import { generateData } from './utils/dataGenerator';

const theme = createTheme({
  palette: {
    mode: 'light',
  },
});

const defaultParams: Parameters = {
  n: 400,
  cmaYield: 0.1,
  gpYield: 0.11,
  wesYield1Tier: 0.37,
  wesYield2Tier: 0.35,
  wesYield3Tier: 0.33,
  expertFee: 165,
  cmaCost: 1529.8,
  gpCost: 1529.8,
  wesCost: 4589.4,
  cmaPPV: 0.7638,
  cmaNPV: 0.9807,
  gpPPV: 0.95,
  gpNPV: 0.99,
  wesPPV: 0.994,
  wesNPV: 0.999,
  aiPrecision: 0.87,
  aiFDR: 0.13,
  aiFOR: 0.05,
  aiNPV: 0.95,
  alpha: 0.5,
  lambda: 0.03
};

function App() {
  const [parameters, setParameters] = useState<Parameters>(defaultParams);
  const [data, setData] = useState(generateData(defaultParams));

  const handleParameterChange = (newParams: Parameters) => {
    setParameters(newParams);
    setData(generateData(newParams));
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
                  onParameterChange={handleParameterChange}
                />
              </Paper>
            </Grid>
            <Grid item xs={12} md={8}>
              <Paper sx={{ p: 2 }}>
                <Plots data={data} />
              </Paper>
            </Grid>
          </Grid>
        </Box>
      </Container>
    </ThemeProvider>
  );
}

export default App; 