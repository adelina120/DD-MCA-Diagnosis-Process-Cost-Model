import React from 'react';
import {
  Box,
  Typography,
  Slider,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  Divider,
  SelectChangeEvent,
  Grid
} from '@mui/material';
import { Parameters } from '../types';

interface ParameterControlsProps {
  parameters: Parameters;
  onParametersChange: (params: Parameters) => void;
  xAxis?: keyof Parameters;
}

const ParameterControls: React.FC<ParameterControlsProps> = ({
  parameters,
  onParametersChange,
  xAxis
}) => {
  const handleChange = (key: keyof Parameters) => (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    const value = parseFloat(event.target.value);
    if (!isNaN(value)) {
      onParametersChange({
        ...parameters,
        [key]: value
      });
    }
  };

  return (
    <Box>
      <Typography variant="h6" gutterBottom>
        Model Parameters
      </Typography>
      <Grid container spacing={2}>
        <Grid item xs={12}>
          <Typography variant="subtitle1" gutterBottom>
            CMA Parameters
          </Typography>
          <Grid container spacing={2}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="CMA Cost"
                type="number"
                value={parameters.cmaCost}
                onChange={handleChange('cmaCost')}
                disabled={xAxis === 'cmaCost'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="CMA Yield"
                type="number"
                value={parameters.cmaYield}
                onChange={handleChange('cmaYield')}
                disabled={xAxis === 'cmaYield'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="CMA PPV"
                type="number"
                value={parameters.cmaPPV}
                onChange={handleChange('cmaPPV')}
                disabled={xAxis === 'cmaPPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="CMA NPV"
                type="number"
                value={parameters.cmaNPV}
                onChange={handleChange('cmaNPV')}
                disabled={xAxis === 'cmaNPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
          </Grid>
        </Grid>

        <Grid item xs={12}>
          <Typography variant="subtitle1" gutterBottom>
            GP Parameters
          </Typography>
          <Grid container spacing={2}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="GP Cost"
                type="number"
                value={parameters.gpCost}
                onChange={handleChange('gpCost')}
                disabled={xAxis === 'gpCost'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="GP Yield"
                type="number"
                value={parameters.gpYield}
                onChange={handleChange('gpYield')}
                disabled={xAxis === 'gpYield'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="GP PPV"
                type="number"
                value={parameters.gpPPV}
                onChange={handleChange('gpPPV')}
                disabled={xAxis === 'gpPPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="GP NPV"
                type="number"
                value={parameters.gpNPV}
                onChange={handleChange('gpNPV')}
                disabled={xAxis === 'gpNPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
          </Grid>
        </Grid>

        <Grid item xs={12}>
          <Typography variant="subtitle1" gutterBottom>
            WES Parameters
          </Typography>
          <Grid container spacing={2}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES Cost"
                type="number"
                value={parameters.wesCost}
                onChange={handleChange('wesCost')}
                disabled={xAxis === 'wesCost'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES Yield"
                type="number"
                value={parameters.wesYield}
                onChange={handleChange('wesYield')}
                disabled={xAxis === 'wesYield'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES Yield 1 Tier"
                type="number"
                value={parameters.wesYield1Tier}
                onChange={handleChange('wesYield1Tier')}
                disabled={xAxis === 'wesYield1Tier'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES Yield 2 Tier"
                type="number"
                value={parameters.wesYield2Tier}
                onChange={handleChange('wesYield2Tier')}
                disabled={xAxis === 'wesYield2Tier'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES Yield 3 Tier"
                type="number"
                value={parameters.wesYield3Tier}
                onChange={handleChange('wesYield3Tier')}
                disabled={xAxis === 'wesYield3Tier'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES PPV"
                type="number"
                value={parameters.wesPPV}
                onChange={handleChange('wesPPV')}
                disabled={xAxis === 'wesPPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="WES NPV"
                type="number"
                value={parameters.wesNPV}
                onChange={handleChange('wesNPV')}
                disabled={xAxis === 'wesNPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
          </Grid>
        </Grid>

        <Grid item xs={12}>
          <Typography variant="subtitle1" gutterBottom>
            Other Parameters
          </Typography>
          <Grid container spacing={2}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Expert Fee"
                type="number"
                value={parameters.expertFee}
                onChange={handleChange('expertFee')}
                disabled={xAxis === 'expertFee'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Utility of True Positive"
                type="number"
                value={parameters.uTP}
                onChange={handleChange('uTP')}
                disabled={xAxis === 'uTP'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Utility of True Negative"
                type="number"
                value={parameters.uTN}
                onChange={handleChange('uTN')}
                disabled={xAxis === 'uTN'}
                InputLabelProps={{ shrink: true }}
              />        
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Utility of False Positive"
                type="number"
                value={parameters.uFP}
                onChange={handleChange('uFP')}
                disabled={xAxis === 'uFP'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Utility of False Negative"
                type="number"
                value={parameters.uFN}
                onChange={handleChange('uFN')}
                disabled={xAxis === 'uFN'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Number of Years"
                type="number"
                value={parameters.numberOfYears}
                onChange={handleChange('numberOfYears')}
                disabled={xAxis === 'numberOfYears'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
          </Grid>
        </Grid>

        <Grid item xs={12}>
          <Typography variant="subtitle1" gutterBottom>
            AI Parameters
          </Typography>
          <Grid container spacing={2}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="AI Precision"
                type="number"
                value={parameters.aiPrecision}
                onChange={handleChange('aiPrecision')}
                disabled={xAxis === 'aiPrecision'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="AI FDR"
                type="number"
                value={parameters.aiFDR}
                onChange={handleChange('aiFDR')}
                disabled={xAxis === 'aiFDR'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="AI FOR"
                type="number"
                value={parameters.aiFOR}
                onChange={handleChange('aiFOR')}
                disabled={xAxis === 'aiFOR'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="AI NPV"
                type="number"
                value={parameters.aiNPV}
                onChange={handleChange('aiNPV')}
                disabled={xAxis === 'aiNPV'}
                InputLabelProps={{ shrink: true }}
              />
            </Grid>
          </Grid>
        </Grid>
      </Grid>
    </Box>
  );
};

export default ParameterControls; 