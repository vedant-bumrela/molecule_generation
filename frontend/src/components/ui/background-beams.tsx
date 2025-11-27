"use client";
import React from "react";
import { motion } from "framer-motion";
import { cn } from "@/lib/utils";

export const BackgroundBeams = ({ className }: { className?: string }) => {
    return (
        <div
            className={cn(
                "absolute h-full w-full inset-0 bg-neutral-950",
                className
            )}
        >
            <motion.div
                initial={{
                    opacity: 0,
                }}
                animate={{
                    opacity: 1,
                    transition: {
                        duration: 2,
                    },
                }}
                className="absolute inset-0 h-full w-full"
            >
                <div className="absolute h-full w-full bg-neutral-950 [mask-image:radial-gradient(ellipse_at_center,transparent_20%,black)]"></div>
                <div className="absolute inset-0 bg-fixed bg-[radial-gradient(#ffffff33_1px,transparent_1px)] [background-size:16px_16px] [mask-image:radial-gradient(ellipse_at_center,black_70%,transparent_100%)] opacity-20"></div>
                <div className="absolute -left-[10%] top-[-10%] h-[500px] w-[500px] rounded-full bg-purple-500/20 blur-[100px] animate-blob mix-blend-screen filter"></div>
                <div className="absolute -right-[10%] top-[-10%] h-[500px] w-[500px] rounded-full bg-cyan-500/20 blur-[100px] animate-blob animation-delay-2000 mix-blend-screen filter"></div>
                <div className="absolute left-[20%] bottom-[-20%] h-[500px] w-[500px] rounded-full bg-blue-500/20 blur-[100px] animate-blob animation-delay-4000 mix-blend-screen filter"></div>
            </motion.div>
        </div>
    );
};
