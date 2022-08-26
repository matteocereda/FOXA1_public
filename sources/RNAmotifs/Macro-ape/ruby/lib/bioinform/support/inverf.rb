module Math
  def self.inverf(x)
    sign = x < 0 ? -1 : 1
    x = x.abs
    a = 8 / (3*Math::PI) * (Math::PI-3) / (4-Math::PI)
    part0 = ( 2/(Math::PI*a) + (Math.log(1-x*x)) / 2 )**2
    part = -2 / (Math::PI * a) - Math.log(1-x*x)/2 + Math.sqrt(-1/a * Math.log(1-x*x) + part0)
    sign * Math.sqrt(part)
  end
  def inverf(x)
    Math.inverf(x)
  end
end