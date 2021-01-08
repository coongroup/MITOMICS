using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ImputeMissingH3kValues
{
    public static class Global
    {
        [ThreadStatic] public static readonly Random Random = new Random(69);
    }
}
