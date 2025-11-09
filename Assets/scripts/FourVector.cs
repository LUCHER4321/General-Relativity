using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FourVector
{
    public double[] components = new double[4];

    public FourVector(double ct, double x, double y, double z)
    {
        components[0] = ct;
        components[1] = x;
        components[2] = y;
        components[3] = z;
    }

    public FourVector(double ct, Vector3 v3, GeneralRelativity config = null)
    {
        components[0] = ct;
        if (config == null || config.currentMetric == MetricType.Minkowski)
        {
            components[1] = v3.x;
            components[2] = v3.y;
            components[3] = v3.z;
        }
        switch (config.currentMetric)
        {
            case MetricType.Schwarzschild:
                components[1] = Math.Sqrt(v3.x * v3.x + v3.y * v3.y + v3.z * v3.z);
                components[2] = Math.Acos(v3.y / components[1]);
                components[3] = Math.Atan2(v3.x, v3.z);
                break;
            case MetricType.Kerr:
                double rho2 = v3.x * v3.x + v3.y * v3.y + v3.z * v3.z;
                double a = config.GetKerrA();
                double a2 = Math.Pow(a, 2.0);
                double sigma = rho2 - a2;
                double sigma2 = Math.Pow(sigma, 2.0);
                components[1] = Math.Sqrt((sigma + Math.Sqrt(sigma2 + 4 * a2 * v3.z * v3.z)) / 2.0);
                if (components[1] > 0)
                {
                    components[2] = Math.Acos(v3.y / components[1]);
                    components[3] = Math.Atan2(v3.x, v3.z);
                }
                break;
        }
    }

    public Vector3 ToVector(GeneralRelativity config = null)
    {
        if (config == null || config.currentMetric == MetricType.Minkowski)
        {
            return new Vector3((float)components[1], (float)components[2], (float)components[3]);
        }
        switch (config.currentMetric)
        {
            case MetricType.Schwarzschild: return (new Vector3((float)Math.Cos(components[3]), 0, (float)Math.Sin(components[3])) * (float)Math.Sin(components[2]) + Vector3.up * (float)Math.Cos(components[2])) * (float)components[1];
            case MetricType.Kerr:
                double r2 = components[1];
                double a = config.GetKerrA();
                double a2 = Math.Pow(a, 2.0);
                double mu = Math.Sqrt(r2 + a2) * Math.Sin(components[2]);
                return new Vector3((float)Math.Cos(components[3]), 0, (float)Math.Sin(components[3])) * (float)mu + Vector3.up * (float)Math.Cos(components[2]) * (float)components[1];
            default: return ToVector();
        }
    }

    //x • y = g_μν × x^μ × y^ν
    public double Dot(SpaceTimeMetric metric, FourVector other)
    {
        double result = 0;
        for (int mu = 0; mu < 4; mu++)
        {
            for (int nu = 0; nu < 4; nu++)
            {
                result += metric.metric[mu, nu] * components[mu] * other.components[nu];
            }
        }
        return result;
    }
}