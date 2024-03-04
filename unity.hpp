#pragma once
#include "mono/mono.h"

namespace unity 
{
	class c_camera;
	class c_shader;
	class c_material;
	class c_renderer;
	class c_time;
	class c_math;
	class c_screen;
	class c_transform;
	class c_asset_bundle;
}

class game_object {
public:
	static auto create_game_object(std::uintptr_t ptr, managed_system::c_string name) -> void {
		MONO_METHOD(internal_create_game_object_fn, "UnityEngine::GameObject.Internal_CreateGameObject()", -1, void(*)(std::uintptr_t, managed_system::c_string));
		return internal_create_game_object_fn(ptr, name);
	}

	static auto add_component(std::uintptr_t ptr, std::uintptr_t type) -> std::uintptr_t
	{
		MONO_METHOD(add_component_fn, "UnityEngine::GameObject.AddComponent()", -1, std::uintptr_t(*)(std::uintptr_t, std::uintptr_t));
		return add_component_fn(ptr, type);
	}
};

class object
{
public:
	static auto dont_destroy_on_load(std::uintptr_t object) -> void
	{
		MONO_METHOD(dont_destroy_on_load_fn, "UnityEngine::Object.DontDestroyOnLoad()", -1, void(*)(std::uintptr_t));
		return dont_destroy_on_load_fn(object);
	}

	static auto destroy(std::uintptr_t object) -> void
	{
		MONO_METHOD(destroy_fn, "UnityEngine::Object.Destroy()", 1, void(*)(std::uintptr_t));
		return destroy_fn(object);
	}
};

#include "c_asset_bundle.hpp"
#include "c_transform.hpp"
#include "c_camera.hpp"
#include "c_shader.hpp"
#include "c_material.hpp"
#include "c_renderer.hpp"
#include "c_math.hpp"
#include "c_time.hpp"
#include "c_screen.hpp"