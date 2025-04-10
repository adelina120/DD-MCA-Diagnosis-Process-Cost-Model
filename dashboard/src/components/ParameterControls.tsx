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
  SelectChangeEvent
} from '@mui/material';
import { Parameters } from '../types';

interface ParameterControlsProps {
  parameters: Parameters;
  onParameterChange: (params: Parameters) => void;
}

const ParameterControls: React.FC<ParameterControlsProps> = ({
  parameters,
  onParameterChange
}) => {
  const handleTextFieldChange = (field: keyof Parameters) => (
    event: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>
  ) => {
    const value = parseFloat(event.target.value);
    onParameterChange({
      ...parameters,
      [field]: value
    });
  };

  const handleSelectChange = (field: keyof Parameters) => (
    event: SelectChangeEvent<number>
  ) => {
    onParameterChange({
      ...parameters,
      [field]: event.target.value
    });
  };

  return (
    <Box>
      <Typography variant="h6" gutterBottom>
        Model Parameters
      </Typography>
      
      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>CMA Parameters</Typography>
        <TextField
          fullWidth
          label="CMA Yield"
          type="number"
          value={parameters.cmaYield}
          onChange={handleTextFieldChange('cmaYield')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="CMA PPV"
          type="number"
          value={parameters.cmaPPV}
          onChange={handleTextFieldChange('cmaPPV')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="CMA NPV"
          type="number"
          value={parameters.cmaNPV}
          onChange={handleTextFieldChange('cmaNPV')}
          margin="normal"
        />
      </Box>

      <Divider sx={{ my: 2 }} />

      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>Gene Panel Parameters</Typography>
        <TextField
          fullWidth
          label="GP Yield"
          type="number"
          value={parameters.gpYield}
          onChange={handleTextFieldChange('gpYield')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="GP PPV"
          type="number"
          value={parameters.gpPPV}
          onChange={handleTextFieldChange('gpPPV')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="GP NPV"
          type="number"
          value={parameters.gpNPV}
          onChange={handleTextFieldChange('gpNPV')}
          margin="normal"
        />
      </Box>

      <Divider sx={{ my: 2 }} />

      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>WES Parameters</Typography>
        <TextField
          fullWidth
          label="WES Yield (1st Tier)"
          type="number"
          value={parameters.wesYield1Tier}
          onChange={handleTextFieldChange('wesYield1Tier')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="WES Yield (2nd Tier)"
          type="number"
          value={parameters.wesYield2Tier}
          onChange={handleTextFieldChange('wesYield2Tier')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="WES Yield (3rd Tier)"
          type="number"
          value={parameters.wesYield3Tier}
          onChange={handleTextFieldChange('wesYield3Tier')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="WES PPV"
          type="number"
          value={parameters.wesPPV}
          onChange={handleTextFieldChange('wesPPV')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="WES NPV"
          type="number"
          value={parameters.wesNPV}
          onChange={handleTextFieldChange('wesNPV')}
          margin="normal"
        />
      </Box>

      <Divider sx={{ my: 2 }} />

      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>Cost Parameters</Typography>
        <TextField
          fullWidth
          label="Expert Fee"
          type="number"
          value={parameters.expertFee}
          onChange={handleTextFieldChange('expertFee')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="CMA Cost"
          type="number"
          value={parameters.cmaCost}
          onChange={handleTextFieldChange('cmaCost')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="GP Cost"
          type="number"
          value={parameters.gpCost}
          onChange={handleTextFieldChange('gpCost')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="WES Cost"
          type="number"
          value={parameters.wesCost}
          onChange={handleTextFieldChange('wesCost')}
          margin="normal"
        />
      </Box>

      <Divider sx={{ my: 2 }} />

      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>Utility Parameters</Typography>
        <TextField
          fullWidth
          label="Alpha"
          type="number"
          value={parameters.alpha}
          onChange={handleTextFieldChange('alpha')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="Lambda"
          type="number"
          value={parameters.lambda}
          onChange={handleTextFieldChange('lambda')}
          margin="normal"
        />
      </Box>

      <Divider sx={{ my: 2 }} />

      <Box sx={{ mb: 2 }}>
        <Typography gutterBottom>AI Parameters</Typography>
        <TextField
          fullWidth
          label="AI Precision"
          type="number"
          value={parameters.aiPrecision}
          onChange={handleTextFieldChange('aiPrecision')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="AI False Discovery Rate (FDR)"
          type="number"
          value={parameters.aiFDR}
          onChange={handleTextFieldChange('aiFDR')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="AI False Omission Rate (FOR)"
          type="number"
          value={parameters.aiFOR}
          onChange={handleTextFieldChange('aiFOR')}
          margin="normal"
        />
        <TextField
          fullWidth
          label="AI Negative Predictive Value (NPV)"
          type="number"
          value={parameters.aiNPV}
          onChange={handleTextFieldChange('aiNPV')}
          margin="normal"
        />
      </Box>
    </Box>
  );
};

export default ParameterControls; 