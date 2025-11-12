using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FourVector
{
    public double[] components = new double[4];

    public FourVector(double[] components)
    {
        this.components = components;
    }

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
                double a = config.GetKerrA().magnitude;
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
                double r2 = Math.Pow(components[1], 2.0);
                Vector3 av3 = config.GetKerrA();
                double a = av3.magnitude;
                double a2 = Math.Pow(a, 2.0);
                double mu = Math.Sqrt(r2 + a2) * Math.Sin(components[2]);
                return (Vector3.right * (float)Math.Cos(components[3]) + Vector3.forward * (float)Math.Sin(components[3])) * (float)mu + Vector3.up * (float)Math.Cos(components[2]) * (float)components[1];
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

    public enum IndexPosition
    {
        up,
        down
    }

    //x_μ = g_μν × x^ν
    //x^μ = g^μν × x_ν
    public FourVector MoveIndex(SpaceTimeMetric metric, IndexPosition pos = IndexPosition.down)
    {
        double[] comps = new double[4];
        double[,] g = new double[4, 4];
        switch (pos)
        {
            case IndexPosition.up:
                g = metric.reverseMetric;
                break;
            case IndexPosition.down:
                g = metric.metric;
                break;
            default: return this;
        }
        for (int mu = 0; mu < 4; mu++)
        {
            comps[mu] = 0;
            for (int nu = 0; nu < 4; nu++)
            {
                comps[mu] += g[mu, nu] * components[nu];
            }
        }
        return new FourVector(comps);
    }

    //u^μ × u_μ = -c²
    public void TimeSpeed(SpaceTimeMetric metric, IndexPosition pos = IndexPosition.up)
    {
        double[,] g = new double[4, 4];
        switch (pos)
        {
            case IndexPosition.up:
                g = metric.reverseMetric;
                break;
            case IndexPosition.down:
                g = metric.metric;
                break;
            default: return;
        }
        double a0 = g[0, 0];
        double b0 = 0;
        double c0 = GeneralRelativity.c2;
        for (int i = 1; i < 4; i++)
        {
            b0 += g[0, i] * components[i];
            for (int j = 1; j < 4; j++)
            {
                c0 += g[i, j] * components[i] * components[j];
            }
        }
        b0 *= 2;
        components[0] = (-b0 + Math.Sqrt(Math.Pow(b0, 2.0) - 4 * a0 * c0)) / (2 * a0);
    }
}