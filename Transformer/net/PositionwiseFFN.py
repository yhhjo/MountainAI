import torch
import torch.nn as nn

class PositionwiseFFN(nn.Module):
    """PositionWiseFFN from Attention is All You Need, dl2.ai, and Cohen et. al (2022).
        Two-layer MLP applied to each index of the input sequence.
        Parameters
        ----------        
        d_model: dimension of latent space in the model.
        d_ffn_hidden: defaults to 2048. Dimension of hidden layer between MLPs.

    """

    def __init__(self, d_model: int, d_ffn_hidden=2048):
        """Default hidden dimension is 2048"""
        super().__init__()
        self._linear1 = nn.Linear(d_model, d_ffn_hidden)
        self._relu = nn.ReLU()
        self._linear2 = nn.Linear(d_ffn_hidden, d_model)


    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Pass x through the PositionwiseFFN block  Input and output have a shape (d_model, d_ffn_hidden).
        Parameters
        ----------
        x: 
            Input tensor of shape (batch_size, input_len, d_model)
        Returns
        ----------
        x:
            Output tensor of shape (batch_size, input_len, d_model)
        """
        x = self._linear1(x)
        x = self._relu(x)
        return self._linear2(x)