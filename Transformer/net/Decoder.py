import torch
import torch.nn as nn
from MultiHeadAttention import MultiHeadAttention
from PositionwiseFFN import PositionwiseFFN
from DotProductAttention import MASKS

class DecoderBlock(nn.Module):
    """ Decoder Block from Attention is All You Need, d2l.ai, and Cohen et. al (2022)
        The goal of the decoder is to predict one target output token at a time given the encoder's output 
        and everything that it has previously predicted. This requires self-attention, encoder-decoder attention,
        and attention masking to preserve autoregression. 

    Layer 1: MHA self attention (with dropout and mask)
            + residual
    Layer 2: Layer norm
    Layer 3: MHA encoder-decoder attention (query = x, k, v = encoder output)  
            + residual
    layer 4: Layer norm
    Layer 5: PositionwiseFFN (with dropout)
            + residual
    Layer 6: Layer norm
    

    Parameters
    ----------
    d_model: 
        Dimension of model's latent space
    q:
        Query and key size
    v:
        Value size
    h: 
        Number of heads
    dropout: 
        Dropout probability for PFFN and MHA blocks
    """
    def __init__(self, d_model, q, v, h, dropout):
        super().__init__()

        self._selfAttention = MultiHeadAttention(d_model, q, v, h, dropout)
        self._layerNorm1 = nn.LayerNorm(d_model)
        self._encoderDecoderAttention = MultiHeadAttention(d_model, q, v, h, dropout)
        self._layerNorm2 = nn.LayerNorm(d_model)
        self._PFFN = PositionwiseFFN(d_model, dropout=dropout)
        self._layerNorm3 = nn.LayerNorm(d_model)

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        """
        Parameters
        ----------
        x: 
            tensor of shape (batch_size, K, d_model)
        context: 
            tensor of shape (batch_size, K, d_model) from the encoding layer output
        Returns
        ----------
        x: 
            tensor of shape (batch_size, K, d_model):
        """
        mask = MASKS.get("future")

        x = self._selfAttention(q=x, k=x, v=x, mask=mask) + x 
        x = self._layerNorm1(x)
        x = self._encoderDecoderAttention(q=x, k=context, v=context) + x
        x = self._layerNorm2(x)
        x = self._PFFN(x) + x
        x = self._layerNorm3(x)
        return x
