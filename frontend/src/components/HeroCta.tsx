import { Button } from "@/components/ui/button";
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card";
import { ArrowRight, FlaskConical, Cpu, Microscope, Leaf } from "lucide-react";
import React from "react";
import { useNavigate } from "react-router-dom";

const FeatureCard: React.FC<{ icon: React.ReactNode; title: string; description: string }> = ({
  icon,
  title,
  description,
}) => (
  <Card className="group relative overflow-hidden border-white/10 bg-white/5 backdrop-blur-lg hover:bg-white/10 transition-all duration-500 hover:scale-[1.02] hover:shadow-2xl hover:shadow-primary/20">
    <div className="absolute inset-0 bg-gradient-to-br from-primary/10 via-transparent to-transparent opacity-0 group-hover:opacity-100 transition-opacity duration-500" />
    <CardHeader className="items-center pb-3 relative z-10">
      <div className="p-4 bg-primary/20 text-primary-foreground rounded-2xl mb-4 group-hover:scale-110 transition-transform duration-500 shadow-glow">
        {icon}
      </div>
      <CardTitle className="text-xl font-bold text-white tracking-tight">{title}</CardTitle>
    </CardHeader>
    <CardContent className="px-6 pb-6 relative z-10">
      <CardDescription className="text-base text-gray-300 text-center leading-relaxed">
        {description}
      </CardDescription>
    </CardContent>
  </Card>
);

export const HeroCta = () => {
  const navigate = useNavigate();

  return (
    <div className="relative overflow-hidden pt-16 pb-20 md:pt-24 md:pb-32">
      {/* FIXED Background Layer - Subtle and Elegant */}
      <div className="fixed inset-0 pointer-events-none z-0">
        {/* Large gradient orbs - REDUCED OPACITY */}
        <div className="absolute top-[10%] left-[5%] w-[700px] h-[700px] bg-gradient-to-br from-primary/35 via-purple-500/25 to-transparent rounded-full blur-[90px] animate-pulse" />
        <div className="absolute bottom-[10%] right-[5%] w-[600px] h-[600px] bg-gradient-to-tl from-secondary/30 via-cyan-500/20 to-transparent rounded-full blur-[80px] animate-pulse" style={{ animationDelay: '1s' }} />
        <div className="absolute top-[30%] right-[20%] w-[500px] h-[500px] bg-gradient-to-br from-violet-500/25 via-pink-500/15 to-transparent rounded-full blur-[70px] animate-pulse" style={{ animationDelay: '2s' }} />

        {/* Glowing particles - SOFTER GLOW */}
        <div className="absolute top-[20%] left-[15%] w-3 h-3 bg-primary/80 rounded-full shadow-[0_0_15px_3px_rgba(124,58,237,0.6)] animate-pulse" style={{ animationDelay: '0.5s' }} />
        <div className="absolute top-[30%] right-[20%] w-2.5 h-2.5 bg-secondary/80 rounded-full shadow-[0_0_12px_3px_rgba(6,182,212,0.6)] animate-pulse" style={{ animationDelay: '1.2s' }} />
        <div className="absolute top-[50%] left-[25%] w-2.5 h-2.5 bg-violet-400/80 rounded-full shadow-[0_0_12px_3px_rgba(167,139,250,0.6)] animate-pulse" style={{ animationDelay: '0.8s' }} />
        <div className="absolute bottom-[25%] right-[30%] w-2.5 h-2.5 bg-cyan-400/80 rounded-full shadow-[0_0_12px_3px_rgba(34,211,238,0.6)] animate-pulse" style={{ animationDelay: '1.5s' }} />
        <div className="absolute top-[60%] right-[45%] w-2.5 h-2.5 bg-pink-400/80 rounded-full shadow-[0_0_12px_3px_rgba(244,114,182,0.6)] animate-pulse" style={{ animationDelay: '2s' }} />
      </div>

      {/* Hero Section */}
      <header className="relative py-24 px-4 sm:px-6 lg:px-8 max-w-7xl mx-auto z-10">
        <div className="relative max-w-4xl mx-auto text-center z-10">
          <div className="inline-flex items-center justify-center gap-3 mb-8 px-6 py-2 rounded-full bg-white/5 border border-white/10 backdrop-blur-md animate-fade-in">
            <FlaskConical className="w-5 h-5 text-primary" />
            <span className="text-sm font-medium text-primary-foreground/90">Next-Gen Molecular Docking</span>
          </div>

          <h1 className="text-6xl sm:text-7xl lg:text-8xl font-extrabold text-white tracking-tight mb-8 leading-tight">
            AI-Enabled <br />
            <span className="text-transparent bg-clip-text bg-gradient-to-r from-primary via-purple-400 to-secondary animate-gradient-x">
              Discovery
            </span>
          </h1>

          <p className="text-xl sm:text-2xl text-gray-300 max-w-3xl mx-auto leading-relaxed mb-12 font-light">
            Accelerate drug discovery with cutting-edge machine learning models and high-throughput computational chemistry.
          </p>

          <div className="flex flex-col sm:flex-row justify-center items-center gap-6">
            <Button
              size="lg"
              className="h-14 px-8 text-lg bg-primary hover:bg-primary/90 text-white shadow-glow hover:shadow-lg hover:shadow-primary/40 transition-all duration-300 rounded-full"
              onClick={() => navigate("/docking")}
            >
              Start Docking Now
              <ArrowRight className="ml-2 h-5 w-5" />
            </Button>
            <Button
              size="lg"
              className="h-14 px-8 text-lg bg-gradient-to-r from-green-500 to-emerald-600 hover:from-green-600 hover:to-emerald-700 text-white shadow-lg hover:shadow-green-500/40 transition-all duration-300 rounded-full"
              onClick={() => navigate("/hdi")}
            >
              <Leaf className="mr-2 h-5 w-5" />
              Herb-Drug Interaction
            </Button>
            <Button
              variant="outline"
              size="lg"
              className="h-14 px-8 text-lg bg-transparent border-white/20 text-white hover:bg-white/10 hover:border-white/40 backdrop-blur-sm transition-all duration-300 rounded-full"
            >
              View Documentation
            </Button>
          </div>
        </div>
      </header>

      {/* Features Section */}
      <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 mt-12 relative z-10">
        <div className="grid grid-cols-1 md:grid-cols-3 gap-8">
          <FeatureCard
            icon={<Cpu className="w-8 h-8" />}
            title="ML Co-Modeling"
            description="Leverage Equibind and GNINA for superior pose prediction and rescoring accuracy."
          />
          <FeatureCard
            icon={<Microscope className="w-8 h-8" />}
            title="High-Throughput"
            description="Process hundreds of molecules efficiently using our optimized batch processing engine."
          />
          <FeatureCard
            icon={<FlaskConical className="w-8 h-8" />}
            title="Advanced Methods"
            description="Support for classical VINA docking alongside our proprietary AI-enhanced protocols."
          />
        </div>
      </div>
    </div>
  );
};
