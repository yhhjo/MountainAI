import torch
import torch.nn as nn
import torch.nn.functional as F
from DotProductAttention import DotProductAttention

class MultiHeadAttention(nn.Module):
    def __init__(self, d_model: int, q: int, v: int, h: int, dropout = 0) -> torch.Tensor:
        """MHA block from Attention is All You Need, d2l.ai, and Cohen et. al (2022). 
        Note: Parallelization is only possible if q = k = v = d_model/h.
        
        Parameters
        ----------
        d_model:
            Dimension of model's latent space
        q:
            Query and key size
        v:
            value size
        h: 
            Number of heads
        """
        super().__init__()

        self._h = h
        self._attention = DotProductAttention(dropout)
        
        # Query, key, value, and output weights
        self._WQ = nn.Linear(d_model, self._h*q)
        self._WK = nn.Linear(d_model, self._h*q)
        self._WV = nn.Linear(d_model, self._h*v)
        self._WO = nn.Linear(self._h*v, d_model)

    def forward(self, q, k, v, mask = None) -> torch.Tensor:
        """Compute query, key, and value matrices to feed into scaled dot product attention. Matrices are divided into size
        (d_model/_h) along the batch dimension to allow parallel computation of multiple attention scores. 
        
        To preserve the autoregressive property, and mask is applied over subsequent tokens (tokens not yet produced). 


        Parameters 
        ----------
        q:
            Query tensor shape (batch_size, K, d_model)
        k:
            Key tensor shape (batch_size, K, d_model)
        v:
            Value tensor shape (batch_size, K, d_model)
        mask:
            None or "forward"

        Returns
        ----------
        Self attention:
            Tensor shape (batch_size, K, d_model)
        """
        # Compute Q, K and V, concatenate heads on batch dimension for parallel processing
        queries = torch.cat(self._WQ(q).chunk(self._h, dim=-1), dim=0)
        keys = torch.cat(self._WK(k).chunk(self._h, dim=-1), dim=0)
        values = torch.cat(self._WV(v).chunk(self._h, dim=-1), dim=0)

        # Save for visualization? 
        attention = self._attention(queries, keys, values, mask)

        # Concatenate all the heads
        attention = torch.cat(attention.chunk(self._h, dim=0), dim=-1)

        # One last linear transformation
        attention =  self._WO(attention)
        return attention
