import React, { useState, useEffect } from 'react';
import axios from 'axios';
import {
  Slider,
  Typography,
  Button,
  LinearProgress,
  Container,
  Box,
} from '@mui/material';

function App() {
  const [omegaValues, setOmegaValues] = useState<number[]>([]);
  const [currentOmega, setCurrentOmega] = useState<number | undefined>(
    undefined
  );
  const [imageSrc, setImageSrc] = useState<string | null>(null);
  const [uploading, setUploading] = useState<boolean>(false);

  const decimalPlaces = 7; // Replace 'n' with the desired number of decimal places

  useEffect(() => {
    if (currentOmega !== undefined) {
      fetchImage(currentOmega);
    }
  }, [currentOmega]);

  const fetchOmegas = async () => {
    try {
      const response = await axios.get('http://127.0.0.1:5000/get_omegas');
      const omegas: number[] = response.data.omegas.map((omega: number) =>
        Number(omega.toFixed(decimalPlaces))
      );
      setOmegaValues(omegas);
      setCurrentOmega(omegas[0]); // Initialize currentOmega with the first omega value
    } catch (error) {
      console.error('Error fetching omega values:', error);
    }
  };

  const fetchImage = async (omega: number) => {
    try {
      const response = await axios.get(
        `http://127.0.0.1:5000/get_image/${omega}`,
        {
          responseType: 'blob',
        }
      );
      const url = URL.createObjectURL(response.data);
      setImageSrc(url);
    } catch (error) {
      console.error('Error fetching image:', error);
    }
  };

  const handleSliderChange = (event: Event, newValue: number | number[]) => {
    const roundedValue = Number((newValue as number).toFixed(decimalPlaces));
    setCurrentOmega(roundedValue);
  };

  const handleFileUpload = async (
    event: React.ChangeEvent<HTMLInputElement>
  ) => {
    if (event.target.files && event.target.files[0]) {
      setUploading(true);
      const file = event.target.files[0];
      const formData = new FormData();
      formData.append('file', file);

      try {
        await axios.post('http://127.0.0.1:5000/upload', formData, {
          headers: {
            'Content-Type': 'multipart/form-data',
          },
        });
        await fetchOmegas();
      } catch (error) {
        console.error('Error uploading file:', error);
      } finally {
        setUploading(false);
      }
    }
  };

  return (
    <Container maxWidth="md" style={{ padding: '2rem' }}>
      <Typography variant="h4" gutterBottom>
        Visualization of Graph network with changing ω
      </Typography>
      <input
        accept=".zip"
        style={{ display: 'none' }}
        id="contained-button-file"
        type="file"
        onChange={handleFileUpload}
      />
      <label htmlFor="contained-button-file">
        <Button
          variant="contained"
          color="primary"
          component="span"
          disabled={uploading}
        >
          Upload ZIP File
        </Button>
      </label>
      {uploading && <LinearProgress style={{ marginTop: '1rem' }} />}
      {omegaValues.length > 0 && currentOmega !== undefined && (
        <Box mt={4}>
          <Slider
            value={currentOmega}
            onChange={handleSliderChange}
            step={null}
            marks={omegaValues.map((omega) => ({
              value: omega,
              label: omega.toFixed(decimalPlaces),
            }))}
            min={Math.min(...omegaValues)}
            max={Math.max(...omegaValues)}
            valueLabelDisplay="auto"
          />
          <Typography variant="h6">
              ω: {currentOmega.toFixed(decimalPlaces)}
          </Typography>
          {imageSrc && (
            <img
              src={imageSrc}
              alt={`Omega ${currentOmega}`}
              style={{ marginTop: '2rem', maxWidth: '100%' }}
            />
          )}
        </Box>
      )}
    </Container>
  );
}

export default App;
