using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;


namespace helperFunctions
{
    public static class HelperFunctions
    {
        public static void Print(string s)
        {
            Console.WriteLine(s);
        }
        public static bool TryParseValue(string value, Type targetType, out object parsedValue)
        {
            if (targetType == typeof(double))
            {
                if (double.TryParse(value, out double doubleValue))
                {
                    parsedValue = doubleValue;
                    return true;
                }
            }
            else if (targetType == typeof(int))
            {
                if (int.TryParse(value, out int intValue))
                {
                    parsedValue = intValue;
                    return true;
                }
            }
            else if (targetType == typeof(bool))
            {
                if (bool.TryParse(value, out bool boolValue))
                {
                    parsedValue = boolValue;
                    return true;
                }
            }
            else if(targetType == typeof(Vector3))
            {
                if (TryParseVector3(value, out Vector3 boolValue))
                {
                    parsedValue = boolValue;
                    return true;
                }
            }
            else if (targetType == typeof(Vector2))
            {
                if (TryParseVector2(value, out Vector2 boolValue))
                {
                    parsedValue = boolValue;
                    return true;
                }
            }
            // Add more handling for other types as needed

            parsedValue = null;
            return false;
        }
        public static string deleteSymbolsWithSpace(string value)
        {
            return value.Replace('<', ' ').Replace('>', ' ').Replace('(', ' ').Replace(')', ' ').Replace('{', ' ').Replace('}', ' ').Replace(" ", "");
        }
        public static bool TryParseVector3(string value, out Vector3 parsedValue)
        {
            parsedValue = new Vector3(0,0,0);
            string s = deleteSymbolsWithSpace(value);
            string[] components = s.Split(',');

            if (components.Length == 3 &&
                float.TryParse(components[0], out float x) &&
                float.TryParse(components[1], out float y) &&
                float.TryParse(components[2], out float z))
            {
                parsedValue = new Vector3(x, y, z);
                return true;
            }

            return false;
        }
        public static bool TryParseVector2(string value, out Vector2 parsedValue)
        {
            parsedValue = new Vector2(0, 0);

            string[] components = value.Split(',');

            if (components.Length == 3 &&
                float.TryParse(components[0], out float x) &&
                float.TryParse(components[1], out float y))
            {
                parsedValue = new Vector2(x, y);
                return true;
            }

            return false;
        }
    }
}
