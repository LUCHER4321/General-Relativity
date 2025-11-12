using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[System.Serializable]
public class SpaceTimeMetric
{
    public double[,] metric = new double[4, 4]; //g_μν
    public double[,] reverseMetric = new double[4, 4]; //g^μν
    public double[,,] christoffelSymbols = new double[4, 4, 4]; //Γ^μ_νλ
    public SpaceTimeMetric()
    {
        //Minkowski's metric
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                metric[i, j] = (i == j) ? (i == 0 ? -1 : 1) : 0;
                reverseMetric[i, j] = metric[i, j];
                for (int k = 0; j < 4; k++)
                {
                    christoffelSymbols[i, j, k] = 0;
                }
            }
        }
    }
}