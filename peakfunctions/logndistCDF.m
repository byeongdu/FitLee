function cdf = logndistCDF(r0, sig, x)
a = r0;
v = sig;
%p = logncdf(x, log(r0), sig);

% Copyright (C) 1995, 1996, 1997  Kurt Hornik
%
% This file is part of Octave.
%
% Octave is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2, or (at your option)
% any later version.
%
% Octave is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, write to the Free
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA.

% -*- texinfo -*-
% @deftypefn {Function File} {} lognormal_cdf (@var{x}, @var{a}, @var{v})
% For each element of @var{x}, compute the cumulative distribution
% function (CDF) at @var{x} of the lognormal distribution with
% parameters @var{a} and @var{v}.  If a random variable follows this
% distribution, its logarithm is normally distributed with mean
% @code{log (@var{a})} and variance @var{v}.
%
% Default values are @var{a} = 1, @var{v} = 1.
% @end deftypefn

% Author: KH <Kurt.Hornik@ci.tuwien.ac.at>
% Description: CDF of the log normal distribution

  if (~((nargin == 1) || (nargin == 3)))
    usage ('lognormal_cdf (x, a, v)');
  end

  if (nargin == 1)
    a = 1;
    v = 1;
  end

  % The following "straightforward" implementation unfortunately does
  % not work (because exp (Inf) -> NaN etc):
  % cdf = normal_cdf (log (x), log (a), v);
  % Hence ...

  if (~isscalar (a) || ~isscalar (v))
    [retval, x, a, v] = common_size (x, a, v);
    if (retval > 0)
      error ('lognormal_cdf: x, a and v must be of common size or scalars');
    end
  end

  cdf = zeros (size (x));

  k = find (isnan (x) | ~(a > 0) | ~(a < Inf) | ~(v > 0) | ~(v < Inf));
  if (any (k))
    cdf(k) = NaN;
  end

  k = find ((x == Inf) & (a > 0) & (a < Inf) & (v > 0) & (v < Inf));
  if (any (k))
    cdf(k) = 1;
  end

  k = find ((x > 0) & (x < Inf) & (a > 0) & (a < Inf) & (v > 0) & (v < Inf));
  if (any (k))
    if (isscalar (a) && isscalar (v))
      cdf(k) = stdnormal_cdf ((log (x(k)) - log (a)) / sqrt (v));
    else
      cdf(k) = stdnormal_cdf ((log (x(k)) - log (a(k))) ./ sqrt (v(k)));
    end
   end
end

function cdf = stdnormal_cdf (x)

  if (nargin ~= 1)
    usage ('stdnormal_cdf (x)');
  end

  sz = size (x);
  if (numel(x) == 0)
    error ('stdnormal_cdf: x must not be empty');
  end

  cdf = (ones (sz) + erf (x / sqrt (2))) / 2;

end

