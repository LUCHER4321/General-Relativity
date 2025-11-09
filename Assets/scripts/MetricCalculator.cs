using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class MetricCalculator : MonoBehaviour
{
    public static SpaceTimeMetric CalculateMetric(FourVector position, GeneralRelativity config)
    {
        switch (config.currentMetric)
        {
            case MetricType.Minkowski: return MinkowskiMetric();
            case MetricType.Schwarzschild: return SchwarzschildMetric(position, config.rb.mass);
            case MetricType.Kerr: return KerrMetric(position, config.rb.mass, config.GetKerrA());
            default: return MinkowskiMetric();
        }
    }
    private static SpaceTimeMetric MinkowskiMetric()
    {
        return new SpaceTimeMetric();
    }

    private static SpaceTimeMetric SchwarzschildMetric(FourVector pos, double mass)
    {
        SpaceTimeMetric metric = new SpaceTimeMetric();
        if (pos.components[1] <= 0)
        {
            return metric;
        }
        double rs = 2.0 * GeneralRelativity.gravitationalConstant * mass / GeneralRelativity.c2;
        double schwarzschildFactor = 1.0 - rs / pos.components[1];
        double r2 = Math.Pow(pos.components[1], 2.0);
        metric.metric[0, 0] = -schwarzschildFactor;
        metric.metric[1, 1] = 1 / schwarzschildFactor;
        metric.metric[2, 2] = r2;
        metric.metric[3, 3] = r2 * Math.Pow(Math.Sin(pos.components[2]), 2.0);
        for (int i = 0; i < 4; i++)
        {
            metric.reverseMetric[i, i] = 1 / metric.metric[i, i];
        }
        //Γ^0_μν
        metric.christoffelSymbols[0, 0, 1] = rs / (2.0 * r2 * schwarzschildFactor);
        metric.christoffelSymbols[0, 1, 0] = metric.christoffelSymbols[0, 0, 1];
        //Γ^1_μν
        metric.christoffelSymbols[1, 0, 0] = rs * schwarzschildFactor / (2.0 * r2);
        metric.christoffelSymbols[1, 1, 1] = -metric.christoffelSymbols[0, 0, 1];
        metric.christoffelSymbols[1, 2, 2] = rs - pos.components[1];
        metric.christoffelSymbols[1, 3, 3] = metric.christoffelSymbols[1, 2, 2] * Math.Pow(Math.Sin(pos.components[2]), 2.0);
        //Γ^2_μν
        metric.christoffelSymbols[2, 1, 2] = 1 / pos.components[1];
        metric.christoffelSymbols[2, 2, 1] = metric.christoffelSymbols[2, 1, 2];
        metric.christoffelSymbols[2, 3, 3] = -Math.Sin(pos.components[2]) * Math.Cos(pos.components[2]);
        //Γ^3_μν
        metric.christoffelSymbols[3, 1, 3] = metric.christoffelSymbols[2, 1, 2];
        metric.christoffelSymbols[3, 2, 3] = 1 / Math.Tan(pos.components[2]);
        metric.christoffelSymbols[3, 3, 1] = metric.christoffelSymbols[2, 1, 2];
        metric.christoffelSymbols[3, 3, 2] = metric.christoffelSymbols[3, 2, 3];
        return metric;
    }
    private static SpaceTimeMetric KerrMetric(FourVector pos, double mass, double a)
    {
        SpaceTimeMetric metric = new SpaceTimeMetric();
        double r2 = Math.Pow(pos.components[1], 2.0);
        double a2 = Math.Pow(a, 2.0);
        double mu = r2 + a2;
        double rs = 2.0 * GeneralRelativity.gravitationalConstant * mass / GeneralRelativity.c2;
        double sin = Math.Sin(pos.components[2]);
        double sin2 = Math.Pow(sin, 2.0);
        double cos2 = 1 - sin2;
        double cos = Math.Sqrt(cos2);
        double upperSigma = r2 + a2 * cos2;
        double upperDelta = mu - rs * pos.components[1];
        metric.metric[0, 0] = -(1 - rs * pos.components[1] / upperSigma);
        metric.metric[0, 3] = -rs * pos.components[1] * a * sin2 / upperSigma;
        metric.metric[1, 1] = upperSigma / upperDelta;
        metric.metric[2, 2] = upperSigma;
        metric.metric[3, 0] = metric.metric[0, 3];
        metric.metric[3, 3] = (mu - a * metric.metric[0, 3]) * sin2;
        double gFactor = metric.metric[0, 0] * metric.metric[3, 3] * Math.Pow(metric.metric[0, 3], 2.0);
        metric.reverseMetric[0, 0] = metric.metric[3, 3] / gFactor;
        metric.reverseMetric[0, 3] = -metric.metric[0, 3] / gFactor;
        metric.reverseMetric[1, 1] = 1 / metric.metric[1, 1];
        metric.reverseMetric[2, 2] = 1 / metric.metric[2, 2];
        metric.reverseMetric[3, 0] = metric.reverseMetric[0, 3];
        metric.reverseMetric[3, 3] = -metric.metric[0, 0] / gFactor;
        double sigma = r2 - a2 * cos2;
        double upperSigma2 = Math.Pow(upperSigma, 2.0);
        double sincos = sin * cos;
        double cot = cos / sin;
        double a4 = Math.Pow(a, 4.0);
        //Γ^0_μν
        metric.christoffelSymbols[0, 0, 1] = rs * mu * sigma / (2.0 * upperSigma2 * upperDelta);
        metric.christoffelSymbols[0, 0, 2] = -rs * pos.components[1] * a2 * sincos / upperSigma2;
        metric.christoffelSymbols[0, 1, 0] = metric.christoffelSymbols[0, 0, 1];
        metric.christoffelSymbols[0, 1, 3] = rs * a * sin2 * (upperSigma * (a2 - r2) - 2 * mu * r2) / (2.0 * upperSigma2 * upperDelta);
        metric.christoffelSymbols[0, 2, 0] = metric.christoffelSymbols[0, 0, 2];
        metric.christoffelSymbols[0, 2, 3] = rs * pos.components[1] * a * a2 * sin2 * sincos / upperSigma2;
        metric.christoffelSymbols[0, 3, 1] = metric.christoffelSymbols[0, 1, 3];
        metric.christoffelSymbols[0, 3, 2] = metric.christoffelSymbols[0, 2, 3];
        //Γ^1_μν
        metric.christoffelSymbols[1, 0, 0] = rs * upperDelta * sigma / (2.0 * upperSigma * upperSigma2);
        metric.christoffelSymbols[1, 0, 3] = metric.christoffelSymbols[1, 0, 0] * a * sin2;
        metric.christoffelSymbols[1, 1, 1] = (2.0 * pos.components[1] * a2 * sin2 - rs * sigma) / (2.0 * upperSigma * upperDelta);
        metric.christoffelSymbols[1, 1, 2] = -a2 * sincos / upperSigma;
        metric.christoffelSymbols[1, 2, 1] = metric.christoffelSymbols[1, 1, 2];
        metric.christoffelSymbols[1, 2, 2] = -pos.components[1] * upperDelta / upperSigma;
        metric.christoffelSymbols[1, 3, 0] = metric.christoffelSymbols[1, 0, 3];
        metric.christoffelSymbols[1, 3, 3] = -metric.christoffelSymbols[1, 0, 3] * sin2 - 2 * pos.components[1] / upperSigma;
        //Γ^2_μν
        metric.christoffelSymbols[2, 0, 0] = metric.christoffelSymbols[1, 1, 2] * rs * pos.components[1] / upperSigma2;
        metric.christoffelSymbols[2, 0, 3] = metric.christoffelSymbols[2, 0, 0] * mu / a;
        metric.christoffelSymbols[2, 1, 1] = metric.christoffelSymbols[1, 1, 2] / upperDelta;
        metric.christoffelSymbols[2, 1, 2] = -pos.components[1] / upperSigma;
        metric.christoffelSymbols[2, 2, 1] = metric.christoffelSymbols[2, 1, 2];
        metric.christoffelSymbols[2, 3, 0] = metric.christoffelSymbols[2, 0, 3];
        metric.christoffelSymbols[2, 3, 3] = -metric.christoffelSymbols[2, 0, 3] * (upperSigma2 / (rs * pos.components[1]) + (1 + upperSigma / mu) * a2 * sin2);
        //Γ^3_μν
        metric.christoffelSymbols[3, 0, 1] = metric.christoffelSymbols[0, 0, 1] * a / mu;
        metric.christoffelSymbols[3, 0, 2] = -rs * pos.components[1] * a * cot / upperSigma2;
        metric.christoffelSymbols[3, 1, 0] = metric.christoffelSymbols[3, 0, 1];
        metric.christoffelSymbols[3, 1, 3] = (rs * (a4 * sin2 * cos2 / (2.0 * upperSigma2 * upperDelta) - r2 * mu * upperSigma) / (2.0 * upperSigma2) + r2) / upperDelta;
        metric.christoffelSymbols[3, 2, 0] = metric.christoffelSymbols[3, 0, 2];
        metric.christoffelSymbols[3, 2, 3] = cot * (1 + rs * pos.components[1] * a2 * sin2 / upperSigma2);
        metric.christoffelSymbols[3, 3, 1] = metric.christoffelSymbols[3, 1, 3];
        metric.christoffelSymbols[3, 3, 2] = metric.christoffelSymbols[3, 2, 3];
        return metric;
    }
}
