import torch
import torch.nn as nn
from MultiHeadAttention import MultiHeadAttention  
from PositionwiseFFN import PositionwiseFFN

class EncoderBlock(nn.Module):
    """Transformer Encoder Block from Attention is All You Need, Cohen et. all, and d2l.ai 
    Layer 1: MultiHeadAttention (no mask, with dropout)
            + Residual
    Layer 2: Layer normalization 
    Layer 3: PositionwiseFFN (with dropout) 
            + Residual connection
    Layer 4: Layer normalization

    Parameters 
    ----------
    d_model:
        Dimension model's latent space
    q:
        Query tensor size
    v:
        Value tensor size
    h: 
        Number of heads 

    attention:
        Attention mode from DotProductAttention.MASKS. Default is none. 

    dropout:
        Dropout rate for dot product attention module.

        Returns
        ----------
        x:
            Tensor shape (batch_size, K, d_model)

    """

    def __init__(self, d_model: int, q: int, v: int, h: int, dropout: float):
        super().__init__()

        self._q = q
        self._v = v
        self._h = h

        self._SelfAttention = MultiHeadAttention(d_model, q, v, h, dropout)
        self._LayerNorm1 = nn.LayerNorm(d_model)
        self._PFFN = PositionwiseFFN(d_model, d_ffn_hidden=1024, dropout=dropout)
        self._LayerNorm2 = nn.LayerNorm(d_model)

        self._dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """ Propogate through Encoder Block. 
        Parameters
        ----------
        x:
            tensor of shape (batch_size, K, d_model)
        Returns
        ----------
        x: 
            tensor of shape (batch_size, K, d_model)
        """

        residual = x
        x = self._SelfAttention(q=x, k=x, v=x)
        x = self._dropout(x) + residual
        x = self._LayerNorm1(x)

        residual = x
        x = self._PFFN(x)
        x = self._dropout(x) + residual
        x = self._LayerNorm2(x)
        return x