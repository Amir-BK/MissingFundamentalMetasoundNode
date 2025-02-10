#include "MetasoundNodeRegistrationMacro.h"
#include "MetasoundParamHelper.h"
#include "MetasoundPrimitives.h"
#include "MetasoundStandardNodesNames.h"
#include "MetasoundStandardNodesCategories.h"
#include "MetasoundFacade.h"
#include "MetasoundVertex.h"
#include "DSP/AlignedBlockBuffer.h"

#define LOCTEXT_NAMESPACE "unDAWMetasounds_MissingFundamentalNode"

namespace unDAWMetasounds_MissingFundamentalNode
{
    METASOUND_PARAM(InputAudio, "Input Audio", "Input audio buffer")
        METASOUND_PARAM(InputCrossoverFrequency, "Crossover Frequency", "Crossover frequency for the filter")
		METASOUND_PARAM(HarmonicGain, "Harmonic Gain", "Gain for the harmonic")
        METASOUND_PARAM(HarmonicDecay, "Harmonic Decay", "Higher harmonics attenuation")
        METASOUND_PARAM(DynamicResponse, "Dynamic Response", "Envelope tracking speed")
		METASOUND_PARAM(Bypass, "Bypass", "Bypass the effect")
        METASOUND_PARAM(OutputAudio, "Output Audio", "Output audio buffer")
}

namespace Metasound
{
    // Parameter names


    class FMissingFundamentalOperator : public TExecutableOperator<FMissingFundamentalOperator>
    {
    public:
        static const FNodeClassMetadata& GetNodeInfo();
        static const FVertexInterface& GetVertexInterface();
        static TUniquePtr<IOperator> CreateOperator(const FBuildOperatorParams& InParams, FBuildResults& OutResults);

        FMissingFundamentalOperator(const FOperatorSettings& InSettings,
            const FAudioBufferReadRef& InAudioInput,
			const FFloatReadRef& InHarmonicGain,
			const FFloatReadRef& InHarmonicDecay,
			const FFloatReadRef& InDynamicResponse,
			const FBoolReadRef& InBypass,
            const FFloatReadRef& InCrossoverFreq);

        virtual void BindInputs(FInputVertexInterfaceData& InOutVertexData) override;
        virtual void BindOutputs(FOutputVertexInterfaceData& InOutVertexData) override;
        virtual void Reset(const IOperator::FResetParams& InParams);
        void Execute();

    private:
        FAudioBufferReadRef AudioInput;
        FAudioBufferWriteRef AudioOutput;
        FFloatReadRef CrossoverFreq;
		FFloatReadRef HarmonicDecay;
		FFloatReadRef DynamicResponse;
		FBoolReadRef Bypass;

        float LpCoeffB0 = 0.0f, LpCoeffB1 = 0.0f, LpCoeffB2 = 0.0f;
        float LpCoeffA1 = 0.0f, LpCoeffA2 = 0.0f;
        float HpCoeffB0 = 0.0f, HpCoeffB1 = 0.0f, HpCoeffB2 = 0.0f;
        float HpCoeffA1 = 0.0f, HpCoeffA2 = 0.0f;

        float FundamentalFreq = 200.0f;
        float EnvelopeRMS = 0.0f;
        float HarmonicBalance[3] = { 0.8f, 0.5f, 0.3f }; // Weighting for 1st, 3rd, 5th harmonics

        FFloatReadRef HarmonicGain; // New member variable

		float HarmonicGainValue = 0.0f; // New member variable

        struct FilterState {
            float z1 = 0.0f;  // First delay element
            float z2 = 0.0f;  // Second delay element
        };

        FilterState LpState;
        FilterState HpState;
        FilterState HpOriginalState;
		FilterState HpHarmonicsState;

        // Filter state variables
        float LpLastInput = 0.0f;
        float LpLastOutput = 0.0f;
        float HpLastInput = 0.0f;
        float HpLastOutput = 0.0f;
        float SampleRate = 48000.0f;

        float LpLastInput2 = 0.0f;
        float LpLastOutput2 = 0.0f;
        float HpLastInput2 = 0.0f;
        float HpLastOutput2 = 0.0f;

        float DcOffset = 0.0f;
        const float DcBlockCoeff = 0.999f;
        bool bFirstExecute = true;

        void UpdateFilterCoefficients(float Crossover);
		void UpdatePsychoacousticModel(float InputSample);
    };

    FMissingFundamentalOperator::FMissingFundamentalOperator(const FOperatorSettings& InSettings,
        const FAudioBufferReadRef& InAudioInput,
		const FFloatReadRef& InHarmonicGain,
		const FFloatReadRef& InHarmonicDecay,
		const FFloatReadRef& InDynamicResponse,
		const FBoolReadRef& InBypass,
        const FFloatReadRef& InCrossoverFreq)
        : AudioInput(InAudioInput)
		, HarmonicGain(InHarmonicGain)
		, HarmonicDecay(InHarmonicDecay)
		, DynamicResponse(InDynamicResponse)
		, Bypass(InBypass)
        , AudioOutput(FAudioBufferWriteRef::CreateNew(InSettings))
        , CrossoverFreq(InCrossoverFreq)
        , SampleRate(InSettings.GetSampleRate())
    {
        UpdateFilterCoefficients(*CrossoverFreq);
    }

    void FMissingFundamentalOperator::BindInputs(FInputVertexInterfaceData& InOutVertexData)
    {
        using namespace unDAWMetasounds_MissingFundamentalNode;

        InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(InputAudio), AudioInput);
		InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(HarmonicGain), HarmonicGain);
		InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(HarmonicDecay), HarmonicDecay);
		InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(DynamicResponse), DynamicResponse);
		InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(Bypass), Bypass);
        InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(InputCrossoverFrequency), CrossoverFreq);
    }

    void FMissingFundamentalOperator::BindOutputs(FOutputVertexInterfaceData& InOutVertexData)
    {
        using namespace unDAWMetasounds_MissingFundamentalNode;
        InOutVertexData.BindReadVertex(METASOUND_GET_PARAM_NAME(OutputAudio), AudioOutput);
    }

    void FMissingFundamentalOperator::Reset(const IOperator::FResetParams& InParams)
    {
        // Reset filter states
        LpState = FilterState{};
        HpState = FilterState{};
        HpOriginalState = FilterState{};
        HpHarmonicsState = FilterState{};
        DcOffset = 0.0f;
        SampleRate = InParams.OperatorSettings.GetSampleRate();
        UpdateFilterCoefficients(*CrossoverFreq);
        AudioOutput->Zero();
    }

    void FMissingFundamentalOperator::UpdatePsychoacousticModel(float InputSample)
    {
        // 1. FUNDAMENTAL FREQUENCY TRACKING
   // Smoothly track crossover frequency changes
        const float RawFrequency = FMath::Clamp(*CrossoverFreq, 20.0f, 500.0f);
        FundamentalFreq = 0.95f * FundamentalFreq + 0.05f * RawFrequency;

        // 2. EQUAL LOUDNESS COMPENSATION
        // ISO 226 curve approximation with parameter-driven adjustments
        const float FrequencyNorm = FundamentalFreq / 1000.0f;
        const float LoudnessCurve = 1.0f - 0.8f * FMath::Exp(-0.45f * FrequencyNorm);
        const float DynamicWeight = FMath::Clamp(*DynamicResponse, 0.01f, 1.0f);
        const float loudnessWeight = LoudnessCurve * DynamicWeight;

        // 3. HARMONIC BALANCE CONTROL
        // Use HarmonicDecay parameter to control harmonic rolloff
        const float baseGain = *HarmonicGain * loudnessWeight;
        const float decayFactor = FMath::Clamp(*HarmonicDecay, 0.0f, 1.0f);

        HarmonicBalance[0] = baseGain;                  // 1st harmonic (full gain)
        HarmonicBalance[1] = baseGain * decayFactor;    // 3rd harmonic 
        HarmonicBalance[2] = baseGain * decayFactor * 0.5f; // 5th harmonic

        // 4. DYNAMIC ENVELOPE FOLLOWER
        // Convert DynamicResponse to time constants
        const float attackTime = FMath::GetMappedRangeValueClamped(
            FVector2f(0.0f, 1.0f), FVector2f(0.001f, 0.1f), *DynamicResponse);
        const float releaseTime = attackTime * 2.0f; // Standard 1:2 ratio

        const float env = FMath::Abs(InputSample);
        const float attackAlpha = 1.0f - FMath::Exp(-1.0f / (attackTime * SampleRate));
        const float releaseAlpha = 1.0f - FMath::Exp(-1.0f / (releaseTime * SampleRate));

        EnvelopeRMS = env > EnvelopeRMS ?
            FMath::Lerp(EnvelopeRMS, env, attackAlpha) :
            FMath::Lerp(EnvelopeRMS, env, releaseAlpha);
    }

    void FMissingFundamentalOperator::Execute()
    {
        const float* InputBuffer = AudioInput->GetData();
        float* OutputBuffer = AudioOutput->GetData();
        const int32 NumSamples = AudioInput->Num();
        const float CurrentCrossover = FMath::Clamp(*CrossoverFreq, 20.0f, SampleRate / 2.1f);

        // Get parameter values once per block
        const float CurrentGain = *HarmonicGain;
        const float CurrentDecay = *HarmonicDecay;
        const float CurrentDynamics = *DynamicResponse;
        const bool bShouldBypass = *Bypass;

        UpdateFilterCoefficients(CurrentCrossover);

        for (int32 i = 0; i < NumSamples; ++i)
        {
            const float input = InputBuffer[i];

            if (bShouldBypass)
            {
                OutputBuffer[i] = input;
                continue;
            }

            // Update psychoacoustic model with current sample
           

            // 1. Remove bass from original signal (high-pass)
            const float highPassedOriginal = HpCoeffB0 * input + HpOriginalState.z1;
            HpOriginalState.z1 = HpCoeffB1 * input - HpCoeffA1 * highPassedOriginal + HpOriginalState.z2;
            HpOriginalState.z2 = HpCoeffB2 * input - HpCoeffA2 * highPassedOriginal;

            // 2. Isolate bass frequencies (low-pass)
            const float lowPassedBass = LpCoeffB0 * input + LpState.z1;
            UpdatePsychoacousticModel(lowPassedBass);
            LpState.z1 = LpCoeffB1 * input - LpCoeffA1 * lowPassedBass + LpState.z2;
            LpState.z2 = LpCoeffB2 * input - LpCoeffA2 * lowPassedBass;

            // 3. Generate harmonics using psychoacoustic model
            const float rectified = FMath::Abs(lowPassedBass);
            const float squared = rectified * FMath::Sign(lowPassedBass); // Asymmetric shaping

            // Generate odd harmonics with phase continuity
            const float phase = 2.0f * UE_PI * FundamentalFreq * i / SampleRate;
            const float harmonic1 = squared * FMath::Sin(phase * 2.0f);
            const float harmonic3 = squared * FMath::Sin(phase * 3.0f);
            const float harmonic5 = squared * FMath::Sin(phase * 5.0f);

            // Apply weighted harmonic balance
            float harmonics = (harmonic1 * HarmonicBalance[0] +
                harmonic3 * HarmonicBalance[1] +
                harmonic5 * HarmonicBalance[2]) / 3.0f;

            // Dynamic compression
            const float compressedHarmonics = harmonics / FMath::Max(0.0001f, EnvelopeRMS);

            // 4. Remove residual lows from harmonics
            const float filteredHarmonics = HpCoeffB0 * compressedHarmonics + HpHarmonicsState.z1;
            HpHarmonicsState.z1 = HpCoeffB1 * compressedHarmonics - HpCoeffA1 * filteredHarmonics + HpHarmonicsState.z2;
            HpHarmonicsState.z2 = HpCoeffB2 * compressedHarmonics - HpCoeffA2 * filteredHarmonics;

            // 5. Combine signals
            OutputBuffer[i] = highPassedOriginal + filteredHarmonics;
        }
    }


    void FMissingFundamentalOperator::UpdateFilterCoefficients(float Crossover)
    {
        const float SafeCrossover = FMath::Clamp(Crossover, 20.0f, SampleRate / 2.1f);
        const float omega = 2.0f * UE_PI * SafeCrossover / SampleRate;
        const float Q = 1.0f / FMath::Sqrt(2.0f); // Butterworth Q

        // Common coefficients for both filters
        const float alpha = FMath::Sin(omega) / (2.0f * Q);
        const float cosw = FMath::Cos(omega);

        // Low-pass coefficients (2nd order Butterworth)
        const float b0LP = (1.0f - cosw) / 2.0f;
        const float b1LP = 1.0f - cosw;
        const float b2LP = b0LP;
        const float a0LP = 1.0f + alpha;

        LpCoeffB0 = b0LP / a0LP;
        LpCoeffB1 = b1LP / a0LP;
        LpCoeffB2 = b2LP / a0LP;
        LpCoeffA1 = -2.0f * cosw / a0LP;
        LpCoeffA2 = (1.0f - alpha) / a0LP;

        // High-pass coefficients (2nd order Butterworth)
        const float b0HP = (1.0f + cosw) / 2.0f;
        const float b1HP = -(1.0f + cosw);
        const float b2HP = b0HP;
        const float a0HP = 1.0f + alpha;

        HpCoeffB0 = b0HP / a0HP;
        HpCoeffB1 = b1HP / a0HP;
        HpCoeffB2 = b2HP / a0HP;
        HpCoeffA1 = -2.0f * cosw / a0HP;
        HpCoeffA2 = (1.0f - alpha) / a0HP;
    }

    const FVertexInterface& FMissingFundamentalOperator::GetVertexInterface()
    {
        using namespace unDAWMetasounds_MissingFundamentalNode;

        static const FVertexInterface Interface(
            FInputVertexInterface(
                TInputDataVertex<FAudioBuffer>(METASOUND_GET_PARAM_NAME_AND_METADATA(InputAudio)),
				TInputDataVertex<float>(METASOUND_GET_PARAM_NAME_AND_METADATA(HarmonicGain), 1.0f),
				TInputDataVertex<float>(METASOUND_GET_PARAM_NAME_AND_METADATA(HarmonicDecay), 0.7f),
				TInputDataVertex<float>(METASOUND_GET_PARAM_NAME_AND_METADATA(DynamicResponse), 0.1f),
				TInputDataVertex<bool>(METASOUND_GET_PARAM_NAME_AND_METADATA(Bypass), false),
                TInputDataVertex<float>(METASOUND_GET_PARAM_NAME_AND_METADATA(InputCrossoverFrequency), 200.0f)
            ),
            FOutputVertexInterface(
                TOutputDataVertex<FAudioBuffer>(METASOUND_GET_PARAM_NAME_AND_METADATA(OutputAudio))
            )
        );
        return Interface;
    }

    const FNodeClassMetadata& FMissingFundamentalOperator::GetNodeInfo()
    {
        auto CreateNodeClassMetadata = []() -> FNodeClassMetadata
            {
                FNodeClassMetadata Metadata;
                Metadata.ClassName = { "unDAW", "MissingFundamental", "Audio" };
                Metadata.MajorVersion = 1;
                Metadata.MinorVersion = 0;
                Metadata.DisplayName = LOCTEXT("MissingFundamentalDisplayName", "Missing Fundamental Generator");
                Metadata.Description = LOCTEXT("MissingFundamentalNodeDescription", "Adds harmonics using missing fundamental effect");
                Metadata.Author = "Amir Ben-Kiki";
				Metadata.PromptIfMissing = INVTEXT("Missing Plugin!");
                Metadata.DefaultInterface = GetVertexInterface();
                Metadata.CategoryHierarchy = { LOCTEXT("AudioEffectsCategory", "Effects") };
                return Metadata;
            };

        static const FNodeClassMetadata Metadata = CreateNodeClassMetadata();
        return Metadata;
    }

    TUniquePtr<IOperator> FMissingFundamentalOperator::CreateOperator(const FBuildOperatorParams& InParams, FBuildResults& OutResults)
    {
        using namespace unDAWMetasounds_MissingFundamentalNode;;

        const FInputVertexInterfaceData& InputData = InParams.InputData;
        return MakeUnique<FMissingFundamentalOperator>(
            InParams.OperatorSettings,
            InputData.GetOrConstructDataReadReference<FAudioBuffer>(METASOUND_GET_PARAM_NAME(InputAudio), InParams.OperatorSettings),
			InputData.GetOrCreateDefaultDataReadReference<float>(METASOUND_GET_PARAM_NAME(HarmonicGain), InParams.OperatorSettings),
			InputData.GetOrCreateDefaultDataReadReference<float>(METASOUND_GET_PARAM_NAME(HarmonicDecay), InParams.OperatorSettings),
			InputData.GetOrCreateDefaultDataReadReference<float>(METASOUND_GET_PARAM_NAME(DynamicResponse), InParams.OperatorSettings),
			InputData.GetOrCreateDefaultDataReadReference<bool>(METASOUND_GET_PARAM_NAME(Bypass), InParams.OperatorSettings),
            InputData.GetOrCreateDefaultDataReadReference<float>(METASOUND_GET_PARAM_NAME(InputCrossoverFrequency), InParams.OperatorSettings));

    }

    class FMissingFundamentalNode : public FNodeFacade
    {
    public:
        FMissingFundamentalNode(const FNodeInitData& InitData)
            : FNodeFacade(InitData.InstanceName, InitData.InstanceID, TFacadeOperatorClass<FMissingFundamentalOperator>())
        {
        }
    };

    METASOUND_REGISTER_NODE(FMissingFundamentalNode)
}

#undef LOCTEXT_NAMESPACE